from mpl_toolkits.basemap import Basemap
import geopandas as gpd
import sqlite3 as sql
import matplotlib.pyplot as plt
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


def get_agent_ids(cur, archetype):
    """ Gets all agentIds from Agententry table for wanted archetype

        agententry table has the following format:
            SimId / AgentId / Kind / Spec /
            Prototype / ParentID / Lifetime / EnterTime

    Parameters
    ----------
    cur: cursor
        sqlite cursor3
    archetype: str
        agent's archetype specification

    Returns
    -------
    id_list: list
        list of all agentId strings
    """
    agents = cur.execute("SELECT agentid FROM agententry WHERE spec "
                         "LIKE '%" + archetype + "%' COLLATE NOCASE"
                         ).fetchall()

    return list(str(agent['agentid']) for agent in agents)


def get_prototype_ids(cur, prototype):
    """ Returns agentid of a prototype

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor
    prototype: str
        name of prototype

    Returns
    -------
    agent_id: list
        list of prototype agent_ids as strings
    """
    ids = cur.execute('SELECT agentid FROM agententry '
                      'WHERE prototype = "' +
                      str(prototype) + '" COLLATE NOCASE').fetchall()

    return list(str(agent['agentid']) for agent in ids)


def get_archetype_position(cur, archetype):
    pos = cur.execute("SELECT agentid, spec, prototype, latitude, longitude "
                      "FROM agentposition WHERE spec LIKE '%" + archetype +
                      "%' COLLATE NOCASE")

    return {str(agent['agentid']): (agent['spec'][10:], agent['prototype'],
                                    (agent['latitude'], agent['longitude']))
            for agent in pos}


def available_archetypes(cur):
	""" Returns agentid of all CYCAMORE archetypes
	excluding ManagerInst, DeployInst, and GrowthRegion.

    Parameters
    ----------
    cur: sqlite cursor
        sqlite cursor

    Returns
    -------
    archetypes: set
        set of archetypes available for plotting
    """
    blacklist = {'ManagerInst',
                 'DeployInst',
                 'GrowthRegion'}
    archetypes = get_archetype_position(cur, 'cycamore')
    archetypes = {v[0] for k, v in archetypes.items()}
    archetypes -= blacklist
    return archetypes


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
    coordinates = [tup[2] for tup in list(prototypes.values())]
    llcrnrlat = min(coordinates[0])
    llcrnrlon = min(coordinates[1])
    urcrnrlat = max(coordinates[0])
    urcrnrlon = max(coordinates[1])
    
    bounds = [llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon]
    return bounds


def plot_agents(cur):
    shapes = ['o', 'v', 's', 'p', '8', 'H']
    fig = plt.figure(figsize=(20, 15))
    bounds = get_bounds(cur)
    sim_map = Basemap(projection='cyl',
                      llcrnrlat=bounds[0],
                      llcrnrlon=bounds[1],
                      urcrnrlat=bounds[2],
                      urcrnrlon=bounds[3])
    sim_map.drawcoastlines()
    sim_map.drawcountries()

    for i, archetype in enumerate(available_archetypes(cur)):
        archetypes = get_archetype_position(cur, archetype)
        for k, v in archetypes.items():
            plt.scatter(v[2][1], v[2][0], marker=shapes[i], label=v[0])
            plt.annotate(v[1], xy=v[2])

    plt.show()
