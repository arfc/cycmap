from collections import defaultdict
from decimal import Decimal
from functools import partial
from mpl_toolkits.basemap import Basemap
import sqlite3 as sql
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np

RAD_TO_DEG = 180 / np.pi
KG_TO_TONS = 1 / 1000
QUANTITY_TO_LINEWIDTH = 1 / 100
BOUND_ADDITION = 0.05
CAPACITY_TO_MARKERSIZE = 0.5


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
                                "LIKE '%" + arch + "%' COLLATE NOCASE")
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


def get_bounds(cur):
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
    prototypes = get_archetype_position(cur, 'Cycamore')
    latitudes = [info[2] for info in list(prototypes.values())]
    longitudes = [info[3] for info in list(prototypes.values())]
    llcrnrlat = min(latitudes)
    llcrnrlon = min(longitudes)
    urcrnrlat = max(latitudes)
    urcrnrlon = max(longitudes)
    extra_lon = (urcrnrlon - llcrnrlon) * BOUND_ADDITION
    extra_lat = (urcrnrlat - llcrnrlat) * BOUND_ADDITION
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
        marker_size[(row['longitude'],
                     row['latitude'])] += row['value'] * CAPACITY_TO_MARKERSIZE
    return marker_size


def list_transactions(cur):
    agententry = {}
    transaction_dict = defaultdict(float)
    query = cur.execute("SELECT endtime FROM FINISH")
    endtime = [row['endtime'] for row in query][0]
    query = cur.execute("SELECT agentid, entertime, lifetime FROM AGENTENTRY")
    for row in query:
        lifetime = row['lifetime']
        if  lifetime == -1:
            lifetime = endtime
        agententry[row['agentid']] = lifetime
    query = cur.execute("SELECT senderid, receiverid, commodity, quantity "
                        "FROM TRANSACTIONS INNER JOIN RESOURCES ON "
                        "TRANSACTIONS.resourceid = RESOURCES.resourceid")
    for row in query:
        lifetime = agententry[row['receiverid']]
        transaction_dict[(row['senderid'],
                          row['receiverid'],
                          row['commodity'])] += row['quantity'] / lifetime
    return transaction_dict


def get_lons_lats_labels(cur, arch, merge=False):
    pos_dict = get_archetype_position(cur, arch)
    lons = [agent[3] for agent in pos_dict.values()]
    lats = [agent[2] for agent in pos_dict.values()]
    labels = [agent for agent in pos_dict.keys()]
    if merge:
        return merge_overlapping_labels(lons, lats, labels)
    return lons, lats, labels


def find_overlap(list_of_sets):
    tracker = defaultdict(list)
    for idx, item in enumerate(list_of_sets):
        tracker[item].append(idx)
    dups = {key: idx for key, idx in tracker.items() if len(idx) > 1}
    return dups


def merge_overlapping_labels(lons, lats, labels):
    dups = find_overlap([(lon, lat) for lon, lat in zip(lons, lats)])
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
    arrows = defaultdict(float)
    for key in transaction_dict.keys():
        sender_id = str(key[0])
        receiver_id = str(key[1])
        commodity = key[2]
        quantity = transaction_dict[key] * KG_TO_TONS
        sender_coord = (positions[sender_id][3],
                        positions[sender_id][2])
        receiver_coord = (positions[receiver_id][3],
                          positions[receiver_id][2])
        arrows[(sender_coord, receiver_coord, commodity)] += quantity
    return arrows


def plot_basemap(cur):
    bounds = get_bounds(cur)
    sim_map = Basemap(projection='cyl',
                      llcrnrlat=bounds[0],
                      llcrnrlon=bounds[1],
                      urcrnrlat=bounds[2],
                      urcrnrlon=bounds[3])
    sim_map.drawcoastlines()
    sim_map.drawcountries()
    sim_map.drawstates()
    sim_map.fillcontinents(color='white', lake_color='aqua', zorder=0)
    sim_map.drawmapboundary(fill_color='lightblue', zorder=-1)
    sim_map.drawparallels(np.arange(10,70,20),labels=[1,1,0,0])
    sim_map.drawmeridians(np.arange(-100,0,20),labels=[0,0,0,1])
    return sim_map


def resize_legend():
    legend = plt.legend(loc=0)
    for handle in legend.legendHandles:
        handle._sizes = [30]
    plt.show()


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


def plot_transaction(cur, sim_map, archs, positions, transaction_dict):
    for arch in archs:
        arrows = transaction_arrows(cur, arch, positions, transaction_dict)
        for key, value in arrows.items():
            point_a = [key[0][0], key[1][0]]
            point_b = [key[0][1], key[1][1]]
            commodity = key[2]
            linewidth = value * QUANTITY_TO_LINEWIDTH
            sim_map.plot(point_a, point_b, linewidth=linewidth)


def main(sqlite_file):
    cur = get_cursor(sqlite_file)
    archs = available_archetypes(cur)
    transaction_dict = list_transactions(cur)
    cycamore_positions = get_archetype_position(cur, 'Cycamore')
    fig = plt.figure(1, figsize=(30, 20))
    sim_map = plot_basemap(cur)
    for i, arch in enumerate(archs):
        if arch == 'Reactor':
            plot_reactors(cur, sim_map)
        else:
            plot_nonreactors(cur, arch, sim_map)
    plot_transaction(cur, sim_map, archs, cycamore_positions, transaction_dict)
    resize_legend()
    