from collections import defaultdict
from mpl_toolkits.basemap import Basemap
import mpld3
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
    """ Returns a cursor to an sqlite file

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
    """ Returns list of all CYCAMORE archetypes  excluding ManagerInst,
    DeployInst, and GrowthRegion.

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------
    archetypes: list
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
    """ Returns a dictionary of agents of an archetype with their prototype,
    specification, and coordinates

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------
    positions: dict
        dictionary with "
        key = agentid, and
        value = list of prototype, spec, latitude, longitude"
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
    """ Returns a list with upper and lower limits of latitude and longitude
    for use with Basemap.

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------
    bounds: list
        list with
        lower left latitude,
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
    """ Returns a dictionary of reactor coordinates and marker size using its
    power capacity

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------
    marker_size: dict
        dictionary with "
        key = (longitude, latitude) in float, and
        value = marker size in float for matplotlib scatter"
    """
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
    """ Returns a dictionary of transactions between agents

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------
    transaction_dict: dict
        dictionary with "
        key = (senderid, receiverid, commodity), and
        value = total transaction quantity over lifetime"
    """
    agententry = {}
    transaction_dict = defaultdict(float)
    query = cur.execute("SELECT endtime FROM FINISH")
    endtime = [row['endtime'] for row in query][0]
    query = cur.execute("SELECT agentid, entertime, lifetime FROM AGENTENTRY")
    for row in query:
        lifetime = row['lifetime']
        if lifetime == -1:
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
    """ Returns longitude, latititude, and agentid as separate lists. Merges
    agentids with same lon, lat into comma sepparated string if merge=True

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor
    arch: str
        Cycamore archetype
    merge: bool
        if True,  agentids with same lon, lat pair are merged
        if False, agentids with same lon, lat pair are not merged

    Returns
    -------
    lons: list
        list of longitudes
    lats: list
        list of latitudes
    labels: list
        list of agentids
    """
    pos_dict = get_archetype_position(cur, arch)
    lons = [agent[3] for agent in pos_dict.values()]
    lats = [agent[2] for agent in pos_dict.values()]
    labels = [agent for agent in pos_dict.keys()]
    if merge:
        lons, lats, labels = merge_overlapping_labels(lons, lats, labels)
    return lons, lats, labels


def find_overlap(list_of_sets):
    """ Returns a dictionary of duplicate items and their index from a list 
    of sets

    Parameters
    ----------
    list_of_sets: list
        list of sets

    Returns
    -------
    dups: dict
        dictionary with "
        key = list of sets, and
        value = index of duplicate set"
    """
    tracker = defaultdict(list)
    for idx, item in enumerate(list_of_sets):
        tracker[item].append(idx)
    dups = {key: idx for key, idx in tracker.items() if len(idx) > 1}
    return dups


def merge_overlapping_labels(lons, lats, labels):
    """ Returns truncated lons, lats, labels lists with duplicates removed and
    labels with same coordinates merged.

    Parameters
    ----------
    lons: list
        list of agent longitudes
    lats: list
        list of agent latitudes
    labels: list
        list of agentids

    Returns
    -------
    lons: list
        list of agent longitudes without duplicates
    lats: list
        list of agent latitudes without duplicates
    labels: list
        list of agentid strings where agentids of the same coordinates
        are merged 
    """
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
    """ Returns a dictionary of transactions between agents for plotting

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor
    arch: str
        Agent archetype
    positions: dict
        dictionary of Cycamore agent positions with "
        key = agentid, and
        value = list of prototype, spec, latitude, longitude"
    transaction_dict: dict
        dictionary of all transactions with "
        key = (senderid, receiverid, commodity), and
        value = total transaction quantity over lifetime"

    Returns
    -------
    arrows: dict
        dictionary of transactions and average quantity moved during lifetime
        with "
        key = tuple of sender coord, receiver coord, and commodity, and
        value = average quantity moved during lifetime in [MTHM]"
    """
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
    """ Returns a matplotlib basemap for the simulation region

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor
    
    Returns
    -------
    sim_map: matplotlib basemap
        matplotlib basemap
    """
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
    sim_map.drawparallels(np.arange(10, 70, 20), labels=[1, 1, 0, 0])
    sim_map.drawmeridians(np.arange(-100, 0, 20), labels=[0, 0, 0, 1])
    return sim_map


def resize_legend(legend):
    """ Resizes scatter plot legends to the same size

    Parameters
    ----------
    legend: matplotlib legend
        matplotlib legend
    
    Returns
    -------
    """
    for handle in legend.legendHandles:
        handle._sizes = [30]


def plot_reactors(cur, basemap):
    lons, lats, labels = get_lons_lats_labels(cur, 'Reactor', True)
    marker_dict = reactor_markers(cur)
    markers = [marker_dict[(lon, lat)] for lon, lat in zip(lons, lats)]
    reactors = basemap.scatter(lons, lats,
                               alpha=0.4,
                               color='grey',
                               label='Reactor',
                               edgecolors='black',
                               s=markers)
    for i, (lon, lat, label) in enumerate(zip(lons, lats, labels)):
        plt.text(lon, lat, label,
                 fontsize=8,
                 verticalalignment='top',
                 horizontalalignment='center')
    return reactors


def plot_nonreactors(cur, arch, basemap):
    colors = ['b', 'g', 'r', 'c', 'm']
    lons, lats, labels = get_lons_lats_labels(cur, arch, True)
    nonreactors = basemap.scatter(lons, lats,
                                  alpha=0.4,
                                  label=str(arch),
                                  color=colors[len(arch) % 5])
    for lon, lat, label in zip(lons, lats, labels):
        plt.text(lon, lat, label,
                 fontsize=8,
                 verticalalignment='center',
                 horizontalalignment='center')
    return nonreactors


def plot_transaction(cur, sim_map, archs, positions, transaction_dict):
    for arch in archs:
        arrows = transaction_arrows(cur, arch, positions, transaction_dict)
        for key, value in arrows.items():
            point_a = [key[0][0], key[1][0]]
            point_b = [key[0][1], key[1][1]]
            commodity = key[2]
            linewidth = value * QUANTITY_TO_LINEWIDTH
            sim_map.plot(point_a, point_b, linewidth=linewidth, zorder=0)


def main(sqlite_file):
    cur = get_cursor(sqlite_file)
    archs = available_archetypes(cur)
    transaction_dict = list_transactions(cur)
    cycamore_positions = get_archetype_position(cur, 'Cycamore')
    fig = plt.figure(1, figsize=(30, 20))
    sim_map = plot_basemap(cur)
    reactors = []
    nonreactors = []
    for i, arch in enumerate(archs):
        if arch == 'Reactor':
            reactors = plot_reactors(cur, sim_map)
        else:
            nonreactors = plot_nonreactors(cur, arch, sim_map)
    # labels 
    plot_transaction(cur, sim_map, archs, cycamore_positions, transaction_dict)
    legend = plt.legend(loc=0)
    resize_legend(legend)
    legend = plt.legend(loc=0)
    
