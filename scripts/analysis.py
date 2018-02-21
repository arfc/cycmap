from collections import defaultdict
from decimal import Decimal
from functools import partial
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point
import geopandas as gpd
import sqlite3 as sql
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np


def get_cursor(file_name):
    """ Connects and returns a cursor to an sqlite output file

    Parameters
    ----------
    file_name: str
        name of the sqlite file

    Returns
    -------
    sqlite cursor3
    """
    con = sql.connect(file_name)
    con.row_factory = sql.Row
    return con.cursor()


def available_archetypes(cur):
    """ Returns list of all CYCAMORE archetypes
    excluding ManagerInst, DeployInst, and GrowthRegion.

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------list
        list of archetypes available for plotting
    """
    blacklist = {'ManagerInst',
                 'DeployInst',
                 'GrowthRegion'}
    query = cur.execute("SELECT spec FROM agentposition WHERE spec"
                        " LIKE '%cycamore%' COLLATE NOCASE")
    archetypes = {agent['spec'][10:] for agent in query} - blacklist
    return list(archetypes)


def get_archetype_position(cur, archetype):
    """ Makes an sqlite query to AgentPosition table to obtain
    agentid, spec, prototype, latitude, and longitude

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------
    positions: dict
        dictionary with key=[agentid],
        and value=[list of prototype, spec, latitude, longitude]
    """
    positions = {}
    if archetype in ['Cycamore', 'CYCAMORE', 'cycamore']:
        archs = available_archetypes(cur)
        for arch in archs:
            query = cur.execute("SELECT agentid, spec, prototype, latitude, "
                                "longitude FROM agentposition WHERE spec "
                                "LIKE '%" + archetype + "%' COLLATE NOCASE")
            position = {str(agent['agentid']): [agent['prototype'],
                                                agent['spec'][10:],
                                                agent['latitude'],
                                                agent['longitude']]
                        for agent in query}
            positions = {**positions, **position}
    else:
        query = cur.execute("SELECT agentid, spec, prototype, latitude, "
                            "longitude FROM agentposition WHERE spec "
                            "LIKE '%" + archetype + "%' COLLATE NOCASE")
        positions = {str(agent['agentid']): [agent['prototype'],
                                             agent['spec'][10:],
                                             agent['latitude'],
                                             agent['longitude']]
                     for agent in query}
    return positions


def archetype_dataframe(cur, archetype):
    """ Returns agentid of all CYCAMORE archetypes
    excluding ManagerInst, DeployInst, and GrowthRegion.

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------
    positions_df: GeoDataFrame
        GeoDataFrame object with agentid, prototype, spec,
        and Point shape
    """
    query = get_archetype_position(cur, archetype)
    positions_df = {}
    positions_df['agentid'] = []
    positions_df['prototype'] = []
    positions_df['spec'] = []
    positions_df['geometry'] = []
    for k, v in query.items():
        positions_df['agentid'].append(k)
        positions_df['prototype'].append(v[0])
        positions_df['spec'].append(v[1])
        positions_df['geometry'].append(Point(v[2], v[3]))
    return gpd.GeoDataFrame(positions_df)


def positions_dictionary(cur):
    archs = {}
    for archetype in available_archetypes(cur):
        archs = {**archs, **get_archetype_position(cur, archetype)}
    return archs


def get_bounds(cur, archs):
    """ Returns the upper and lower limits of latitude
    and longitude for use with Basemap.

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------
    bounds: list
        list with lower left latitude,
        lower left longitude,
        upper right latitude,
        upper right longitude in decimal degrees
    """
    prototypes = {}
    for arch in archs:
        arch_dict = get_archetype_position(cur, arch)
        prototypes = {**prototypes, **arch_dict}
    latitudes = [info[2] for info in list(prototypes.values())]
    longitudes = [info[3] for info in list(prototypes.values())]
    llcrnrlat = min(latitudes)
    llcrnrlon = min(longitudes)
    urcrnrlat = max(latitudes)
    urcrnrlon = max(longitudes)
    extra_lon = (urcrnrlon - llcrnrlon) * 0.1 / 2
    extra_lat = (urcrnrlat - llcrnrlat) * 0.1 / 2
    llcrnrlat -= extra_lat
    llcrnrlon -= extra_lon
    urcrnrlat += extra_lat
    urcrnrlon += extra_lon
    bounds = [llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon]
    return bounds


def reactor_markers(cur):
    query = cur.execute("SELECT DISTINCT longitude, latitude, value FROM "
                        "TIMESERIESPOWER INNER JOIN AGENTPOSITION ON "
                        "TIMESERIESPOWER.agentid = AGENTPOSITION.agentid "
                        "EXCEPT "
                        "SELECT DISTINCT longitude, latitude, value FROM "
                        "TIMESERIESPOWER INNER JOIN AGENTPOSITION ON "
                        "TIMESERIESPOWER.agentid = AGENTPOSITION.agentid "
                        "WHERE value = 0")
    marker_size = defaultdict(float)
    for row in query:
        marker_size[(row['longitude'], row['latitude'])] += row['value'] / 2
    return marker_size


def list_transactions(cur):
    agententry = {}
    query = cur.execute("SELECT endtime FROM FINISH")
    endtime = [row['endtime'] for row in query][0]
    query = cur.execute("SELECT agentid, entertime, lifetime FROM AGENTENTRY")
    for row in query:
        agentid = row['agentid']
        lifetime = row['lifetime']
        if lifetime == -1:
            lifetime = endtime
        agententry[agentid] = lifetime
    query = cur.execute("SELECT senderid, receiverid, commodity, quantity "
                        "FROM TRANSACTIONS INNER JOIN RESOURCES ON "
                        "TRANSACTIONS.resourceid = RESOURCES.resourceid")
    transaction_dict = defaultdict(float)
    for row in query:
        senderid = row['senderid']
        receiverid = row['receiverid']
        quantity = row['quantity']
        commodity = row['commodity']
        lifetime = agententry[receiverid]
        transaction_dict[(senderid, receiverid, commodity)
                         ] += quantity / (lifetime * 1000)
    return transaction_dict


def get_lons_lats_labels(cur, arch, merge=False):
    pos_dict = get_archetype_position(cur, arch)
    lons = [agent[3] for agent in pos_dict.values()]
    lats = [agent[2] for agent in pos_dict.values()]
    labels = [agent for agent in pos_dict.keys()]
    if merge:
        return merge_overlapping_labels(lons, lats, labels)
    return lons, lats, labels


def find_overlap(lons, lats, labels):
    coords = [(lon, lat) for lon, lat in zip(lons, lats)]
    tracker = defaultdict(list)
    for idx, item in enumerate(coords):
        tracker[item].append(idx)
    dups = {key: idx for key, idx in tracker.items() if len(idx) > 1}
    return dups, lons, lats, labels


def merge_overlapping_labels(lons, lats, labels):
    dups, lons, lats, labels = find_overlap(lons, lats, labels)
    new_label = []
    dup_idxs = sum(dups.values(), [])
    keep_one_idx_list = [idx[0] for idx in dups.values()]
    for idx in keep_one_idx_list:
        dup_idxs.remove(idx)
    for idx in sorted(dup_idxs, reverse=True):
        del lons[idx]
        del lats[idx]
    for i, item in enumerate(labels):
        if i in keep_one_idx_list:
            new_label.append(', '.join(map(str, [labels[k]
                                                 for j in dups.values()
                                                 if j[0] == i
                                                 for k in j])))
        elif i in dup_idxs:
            continue
        else:
            new_label.append(item)
    return lons, lats, new_label


def transaction_arrows(cur, arch, positions, transaction_dict):
    lons, lats, labels = get_lons_lats_labels(cur, arch, True)
    rad_to_deg = 180 / np.pi
    arrows = defaultdict(float)
    for key in transaction_dict.keys():
        sender_id = str(key[0])
        receiver_id = str(key[1])
        commodity = str(key[0]) + ' -> ' + str(key[1])
        quantity = transaction_dict[key]
        sender_coord = (positions[sender_id][3],
                        positions[sender_id][2])
        receiver_coord = (positions[receiver_id][3],
                          positions[receiver_id][2])
        theta = np.arctan2(sender_coord[1] - receiver_coord[1],
                           sender_coord[0] - receiver_coord[0]) * rad_to_deg
        trans = ((sender_coord[0] - receiver_coord[0]) * 0.15,
                 (sender_coord[1] - receiver_coord[1]) * 0.15)
        arrow_coord = (sender_coord[0] - trans[0],
                       sender_coord[1] - trans[1])
        arrows[(arrow_coord, theta, commodity)] += quantity
    return arrows


def plot_basemap(cur, fig, archs):
    bounds = get_bounds(cur, archs)
    sim_map = Basemap(projection='cyl',
                      llcrnrlat=bounds[0],
                      llcrnrlon=bounds[1],
                      urcrnrlat=bounds[2],
                      urcrnrlon=bounds[3])
    sim_map.drawcoastlines()
    sim_map.drawcountries()
    sim_map.fillcontinents(color='white', lake_color='aqua', zorder=0)
    sim_map.drawmapboundary(fill_color='lightblue', zorder=-1)
    return sim_map


def plot_reactors(cur, basemap):
    lons, lats, labels = get_lons_lats_labels(cur, 'Reactor', True)
    marker_dict = reactor_markers(cur)
    for i, (lon, lat, label) in enumerate(zip(lons, lats, labels)):
        if i == 0:
            basemap.scatter(lon, lat,
                            alpha=0.4,
                            color='grey',
                            label='Reactor',
                            edgecolors='black',
                            s=marker_dict[(lon, lat)])
        else:
            basemap.scatter(lon, lat,
                            alpha=0.4,
                            color='grey',
                            label='_nolegend_',
                            edgecolors='black',
                            s=marker_dict[(lon, lat)])
        plt.text(lon, lat, label,
                 fontsize=8,
                 verticalalignment='top',
                 horizontalalignment='center')


def plot_nonreactors(cur, arch, basemap):
    colors = ['b', 'g', 'r', 'c', 'm']
    lons, lats, labels = get_lons_lats_labels(cur, arch, True)
    basemap.scatter(lons, lats,
                    alpha=0.4,
                    label=str(arch),
                    color=colors[len(arch) % 5])
    for lon, lat, label in zip(lons, lats, labels):
        plt.text(lon, lat, label,
                 fontsize=8,
                 verticalalignment='center',
                 horizontalalignment='center')


def plot_transaction(cur, archs, positions, transaction_dict):
    textbox_arrow_property = dict(boxstyle='larrow, pad=0.3',
                                  fc='cyan', alpha=0.1,
                                  ec='b',
                                  lw=2)
    for arch in archs:
        arrows = transaction_arrows(cur, arch, positions, transaction_dict)
        for key, value in arrows.items():
            coord = key[0]
            theta = key[1]
            commodity = key[2]
            quantity = value
            label = commodity + ':\n' + \
                "{:.3E}".format(Decimal(quantity)) + ' MTHM/month'
            plt.text(coord[0], coord[1],
                     label,
                     rotation=theta, size=7,
                     ha='center', va='center',
                     bbox=textbox_arrow_property)


def main(sqlite_file):
    cur = get_cursor(sqlite_file)
    # Commonly used items
    archs = available_archetypes(cur)
    transaction_dict = list_transactions(cur)
    cycamore_positions = get_archetype_position(cur, 'Cycamore')
    fig = plt.figure(1, figsize=(30, 20))
    sim_map = plot_basemap(cur, fig, archs)
    for i, arch in enumerate(archs):
        if arch == 'Reactor':
            plot_reactors(cur, sim_map)
        else:
            plot_nonreactors(cur, arch, sim_map)
    plot_transaction(cur, arch, cycamore_positions, transaction_dict)
    legend = plt.legend(loc=0)
    for handle in legend.legendHandles:
        handle._sizes = [30]
    plt.show()
