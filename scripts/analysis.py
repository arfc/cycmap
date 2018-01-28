from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point
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
    archetypes: set
        set of archetypes available for plotting
    """
    query = get_archetype_position(cur, archetype)
    positions_df = {}
    for k, v in query.items():
        positions_df['agentid'] = 

    positions_df = {'agentid': [agent['agentid'] for agent in query],
                    'prototype': [agent['prototype'] for agent in query],
                    'spec': [agent['spec'] for agent in query],
                    'geometry': [Point(agent['latitude'],
                                       agent['longitude']) for agent in query]}
    return gpd.GeoDataFrame(positions_df)


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
    bounds = [llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon]
    return bounds


def plot_agents(cur):
    fig = plt.figure(figsize=(20, 15))
    bounds = get_bounds(cur)
    sim_map = Basemap(projection='cyl',
                      llcrnrlat=bounds[0],
                      llcrnrlon=bounds[1],
                      urcrnrlat=bounds[2],
                      urcrnrlon=bounds[3])
    sim_map.drawcoastlines()
    sim_map.drawcountries()
    archs = {}
    for archetype in available_archetypes(cur):
        archs[archetype] = archetype_dataframe(cur, archetype)
        archs[archetype].plot()
    plt.show()
