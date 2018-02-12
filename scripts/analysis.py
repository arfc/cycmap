from collections import defaultdict
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
    query = cur.execute("SELECT agentid, spec, prototype, latitude, longitude "
                        "FROM agentposition WHERE spec LIKE '%" + archetype +
                        "%' COLLATE NOCASE")
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


def capacity_marker(cur):
    query = cur.execute("SELECT ")


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
    archs = available_archetypes(cur)
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
    power_dict = defaultdict(float)
    for row in query:
        power_dict[(row['longitude'], row['latitude'])] += row['value'] / 5
    return power_dict


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


def get_lons_lats_labels(cur, arch):
    pos_dict = get_archetype_position(cur, arch)
    lons = [agent[3] for agent in pos_dict.values()]
    lats = [agent[2] for agent in pos_dict.values()]
    labels = [agent for agent in pos_dict.keys()]
    return lons, lats, labels


def list_transactions(cur):
    agententry = {}
    query = cur.execute("SELECT agentid, entertime, lifetime FROM AGENTENTRY")
    for row in query:
        agentid = row['agentid']
        agententry[agentid] = lifetime
    query = cur.execute("SELECT senderid, receiverid, commodity, quantity "
                        "FROM TRANSACTIONS INNER JOIN RESOURCES ON "
                        "TRANSACTIONS.resourceid = RESOURCES.resourceid")
    transaction_dict = defaultdict(float)
    for row in query:
        senderid = row['senderid']
        receiverid = row['receiverid']
        commodity = row['commodity']
        lifetime = agententry[receiverid]
        transaction_dict[(senderid, receiverid, commodity)
                         ] += quantity / lifetime
    return transaction_dict


def transaction_arrows(cur, transaction_dict):
    


def plot_agents(cur):
    fig = plt.figure(1, figsize=(30, 20))
    bounds = get_bounds(cur)
    sim_map = Basemap(projection='cyl',
                      llcrnrlat=bounds[0],
                      llcrnrlon=bounds[1],
                      urcrnrlat=bounds[2],
                      urcrnrlon=bounds[3])
    sim_map.drawcoastlines()
    sim_map.drawcountries()
    sim_map.fillcontinents(color='white', lake_color='aqua', zorder=0)
    sim_map.drawmapboundary(fill_color='lightblue', zorder=-1)
    for i, arch in enumerate(available_archetypes(cur)):
        if arch == 'Reactor':
            lons, lats, labels = get_lons_lats_labels(cur, arch)
            lons, lats, labels = merge_overlapping_labels(lons, lats, labels)
            marker_dict = reactor_markers(cur)
            for lon, lat, label in zip(lons, lats, labels):
                plt.text(lon, lat, label,
                         fontsize=8,
                         verticalalignment='top',
                         horizontalalignment='center')
                sim_map.scatter(lon, lat,
                                alpha=0.4,
                                color='grey',
                                edgecolors='black',
                                s=marker_dict[(lon, lat)])
        else:
            lons, lats, labels = get_lons_lats_labels(cur, arch)
            lons, lats, labels = merge_overlapping_labels(lons, lats, labels)
            sim_map.scatter(lons, lats)
            for lon, lat, label in zip(lons, lats, labels):
                plt.text(lon, lat, label,
                         fontsize=8,
                         verticalalignment='center',
                         horizontalalignment='center')
    plt.show()
