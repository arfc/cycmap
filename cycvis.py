from collections import defaultdict
from mpl_toolkits.basemap import Basemap
import sqlite3 as sql
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import sys
import warnings
warnings.filterwarnings('ignore')


class Cycvis():
    rad_to_deg = 180 / np.pi
    kg_to_tons = 1 / 1000
    quantity_to_linewidth = 0.1
    bound_addition = 0.05
    capacity_to_markersize = 0.5

    # default settings for mpl plotting
    figsize = (9.8, 5)
    main_plot_axis_position = [0.05, 0.05, 0.6, 0.9]
    annot_property = {'xy': (0, 0),
                      'xytext': (1.02, 0.9785),
                      'textcoords': 'axes fraction',
                      'bbox': dict(boxstyle="round",
                                   alpha=(0.4),
                                   fc="w")}
    label_property = {'fontsize': 8,
                      'vert_align': 'center',
                      'horz_align': 'center'}

    def __init__(self, file_name):
        self.cur = self.get_cursor(file_name)
        self.archs = self.available_archetypes()
        self.transactions = self.list_transactions()
        self.agent_info = self.get_agent_info()
        self.bounds = self.get_bounds()
        self.fig, self.ax, self.basemap = self.plot_basemap()
        self.reactor_markers = self.reactor_markers()
        self.colors = cm.rainbow(np.linspace(0, 1, len(self.archs)))
        self.mpl_collections = self.plot_archetypes()
        self.init_yr, self.timestep = self.sim_info()

    def get_cursor(self, file_name):
        """ Returns a cursor to an sqlite file

        Parameters
        ----------
        file_name: str
                name of the sqlite file

        Returns
        -------
        cur: sqlite cursor
                sqlite cursor
        """
        con = sql.connect(file_name)
        con.row_factory = sql.Row
        cur = con.cursor()
        return cur

    def available_archetypes(self):
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
        query = ("SELECT spec FROM agentposition"
                 " WHERE spec LIKE '%cycamore%'"
                 " COLLATE NOCASE")
        results = self.cur.execute(query)
        archetypes = {agent['spec'][10:] for agent in results} - blacklist
        return list(archetypes)

    def get_agent_info(self):
        """ Returns a dictionary of agents of an archetype with their prototype,
        specification, and coordinates

        Parameters
        ----------
        cur: sqlite cursor
                sqlite cursor

        Returns
        -------
        agent_info: dict
                dictionary with "
                key = agentid, and
                value = list of prototype, spec, latitude, longitude"
        """
        query = ("SELECT agentposition.agentid, agentposition.spec,"
                 " agentposition.prototype, latitude, longitude,"
                 " entertime, lifetime"
                 " FROM agentposition"
                 " INNER JOIN agententry"
                 " ON agentposition.agentid = agententry.agentid"
                 " WHERE agentposition.spec LIKE '%Cycamore%'"
                 " AND agentposition.spec NOT LIKE '%ManagerInst%'"
                 " AND agentposition.spec NOT LIKE '%DeployInst%'"
                 " AND agentposition.spec NOT LIKE '%GrowthRegion%'"
                 " COLLATE NOCASE")
        results = self.cur.execute(query)
        agent_info = {str(agent['agentid']): [agent['prototype'],
                                              agent['spec'][10:],
                                              agent['latitude'],
                                              agent['longitude'],
                                              agent['entertime'],
                                              agent['lifetime']
                                              ]
                      for agent in results}
        return agent_info

    def get_archetype_info(self, archetype):
        archetype_info = {}
        for k, v in self.agent_info.items():
            spec = v[1]
            if spec == archetype:
                archetype_info[k] = v
        return archetype_info

    def get_bounds(self):
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
        latitudes = [info[2] for info in list(self.agent_info.values())]
        longitudes = [info[3] for info in list(self.agent_info.values())]
        llcrnrlat = min(latitudes)
        llcrnrlon = min(longitudes)
        urcrnrlat = max(latitudes)
        urcrnrlon = max(longitudes)
        extra_lon = (urcrnrlon - llcrnrlon) * self.bound_addition
        extra_lat = (urcrnrlat - llcrnrlat) * 3 * self.bound_addition
        llcrnrlat -= extra_lat
        llcrnrlon -= extra_lon
        urcrnrlat += extra_lat
        urcrnrlon += extra_lon
        bounds = [llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon]
        return bounds

    def reactor_markers(self):
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
        query = ("SELECT DISTINCT longitude, latitude, value"
                 " FROM TIMESERIESPOWER INNER JOIN AGENTPOSITION"
                 " ON TIMESERIESPOWER.agentid = AGENTPOSITION.agentid"
                 " WHERE value != 0")
        results = self.cur.execute(query)
        marker_size = defaultdict(float)
        for row in results:
            marker_size[(row['longitude'],
                         row['latitude'])] += (self.capacity_to_markersize *
                                               row['value'])
        return marker_size

    def sim_info(self):
        """ Returns simulation start year, month, duration and
        timesteps (in numpy linspace).

        Parameters
        ----------
        cur: sqlite cursor
            sqlite cursor

        Returns
        -------
        init_year: int
            start year of simulation
        init_month: int
            start month of simulation
        duration: int
            duration of simulation
        timestep: list
            linspace up to duration
        """
        query = 'SELECT initialyear, initialmonth, duration FROM info'
        results = self.cur.execute(query).fetchone()
        init_year = results['initialyear']
        init_month = results['initialmonth']
        duration = results['duration']
        timestep = np.linspace(init_month, init_month + duration - 1, duration)
        return init_year, timestep

    def list_transactions(self):
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
        query = self.cur.execute("SELECT endtime FROM FINISH")
        endtime = [row['endtime'] for row in query][0]
        query = self.cur.execute(
            "SELECT agentid, entertime, lifetime FROM AGENTENTRY")
        for row in query:
            lifetime = row['lifetime']
            if lifetime == -1:
                lifetime = endtime
            agententry[row['agentid']] = lifetime
        query = ("SELECT senderid, receiverid, commodity, quantity"
                 " FROM TRANSACTIONS INNER JOIN RESOURCES ON"
                 " TRANSACTIONS.resourceid = RESOURCES.resourceid")
        results = self.cur.execute(query)
        for row in results:
            lifetime = agententry[row['receiverid']]
            transaction_dict[(row['senderid'],
                              row['receiverid'],
                              row['commodity'])] += row['quantity'] / lifetime
        return transaction_dict

    def get_lons_lats_labels(self, arch, merge=False):
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
        arch_info = self.get_archetype_info(arch)
        lons = [agent[3] for agent in arch_info.values()]
        lats = [agent[2] for agent in arch_info.values()]
        labels = [agent for agent in arch_info.keys()]
        if merge:
            lons, lats, labels = self.merge_overlapping_labels(lons,
                                                               lats,
                                                               labels)
        return lons, lats, labels

    def find_overlap(self, list_of_sets):
        """ Returns a dictionary of duplicate items and their index from a list
        of sets

        Parameters
        ----------
        list_of_tuples: list
                list of tuples

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

    def merge_overlapping_labels(self, lons, lats, labels):
        """ Returns truncated lons, lats, labels lists with duplicates removed and
        labels with same coordinates merged.

        Parameters
        ----------
        lons: list
                list of agent longitudes
        lats: list
                list of agent latitudes
        labels: list
                list of agent prototype names

        Returns
        -------
        lons: list
                list of agent longitudes without duplicates
        lats: list
                list of agent latitudes without duplicates
        new_label: list
                list of agentid (prototype, spec) where agentids of the same
                coordinates are merged
        """
        dups = self.find_overlap([(lon, lat) for lon, lat in zip(lons, lats)])
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

    def transaction_arrows(self, arch):
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
                dictionary of transactions and average quantity moved during
                lifetime with "
                key = tuple of sender coord, receiver coord, and commodity, and
                value = average quantity moved during lifetime in [MTHM]"
        """
        lons, lats, labels = self.get_lons_lats_labels(arch, True)
        arrows = defaultdict(float)
        for key in self.transactions.keys():
            sender_id = str(key[0])
            receiver_id = str(key[1])
            commodity = key[2]
            quantity = self.transactions[key] * self.kg_to_tons
            sender_coord = (self.agent_info[sender_id][3],
                            self.agent_info[sender_id][2])
            receiver_coord = (self.agent_info[receiver_id][3],
                              self.agent_info[receiver_id][2])
            arrows[(sender_coord, receiver_coord, commodity)] += quantity
        return arrows

    def plot_basemap(self):
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
        fig = plt.figure(figsize=self.figsize)
        ax_main = fig.add_axes(self.main_plot_axis_position)
        basemap = Basemap(ax=ax_main,
                          projection='cyl',
                          llcrnrlat=self.bounds[0],
                          llcrnrlon=self.bounds[1],
                          urcrnrlat=self.bounds[2],
                          urcrnrlon=self.bounds[3],
                          fix_aspect=False,
                          anchor='NW')
        basemap.drawcoastlines(zorder=-15)
        basemap.drawmapboundary(fill_color='lightblue', zorder=-10)
        basemap.fillcontinents(color='white', lake_color='aqua', zorder=-5)
        basemap.drawcountries(zorder=0)
        basemap.drawstates(zorder=0)
        return fig, ax_main, basemap

    def resize_legend(self):
        """ Resizes scatter plot legends to the same size

        Parameters
        ----------
        legend: matplotlib legend
                matplotlib legend

        Returns
        -------
        """
        legend = self.ax.legend(loc='best')
        for handle in legend.legendHandles:
            handle._sizes = [30]

    def plot_reactors(self):
        """ Scatter plot of reactors with reactor capacity as marker size

        Parameters
        ----------
        cur: sqlite cursor
                sqlite cursor
        basemap: matplotlib basemap
                matplotlib basemap

        Returns
        -------
        """
        mpl = {}
        lons, lats, labels = self.get_lons_lats_labels('Reactor', True)
        for i, (lon, lat, label) in enumerate(zip(lons, lats, labels)):
            agents = set(label.split(', '))
            mpl[self.basemap.scatter(lon, lat,
                                     alpha=0.4,
                                     color='grey',
                                     label='Reactor' if i == 0 else '',
                                     edgecolors='black',
                                     s=self.reactor_markers[(lon, lat)],
                                     zorder=5)] = agents
            plt.text(lon, lat, label,
                     fontsize=self.label_property['fontsize'],
                     verticalalignment=self.label_property['vert_align'],
                     horizontalalignment=self.label_property['horz_align'])
        return mpl

    def plot_nonreactors(self, arch, i):
        """ Scatter plot of nonreactors

        Parameters
        ----------
        cur: sqlite cursor
                sqlite cursor
        arch: str
                cycamore archetype
        basemap: matplotlib basemap
                matplotlib basemap

        Returns
        -------
        """
        mpl = {}
        lons, lats, labels = self.get_lons_lats_labels(arch, True)
        for j, (lon, lat, label) in enumerate(zip(lons, lats, labels)):
            agents = set(label.split(', '))
            mpl[self.basemap.scatter(lon, lat,
                                     alpha=0.4, s=200,
                                     label=str(arch) if j == 0 else '',
                                     color=self.colors[i],
                                     zorder=5)] = agents
            plt.text(lon, lat, label,
                     fontsize=self.label_property['fontsize'],
                     verticalalignment=self.label_property['vert_align'],
                     horizontalalignment=self.label_property['horz_align'])
        return mpl

    def plot_transactions(self):
        """ Line plot of transactions during simulation with logarithmic linewidth

        Parameters
        ----------
        cur: sqlite cursor
                sqlite cursor
        fig: matplotlib figure
                matplotlib figure
        archs: list
                list of cycamore archetypes

        Returns
        -------
        """
        for arch in self.archs:
            arrows = self.transaction_arrows(arch)
            for key, value in arrows.items():
                point_a = [key[0][0], key[1][0]]
                point_b = [key[0][1], key[1][1]]
                linewidth = np.log(value * self.quantity_to_linewidth)
                self.ax.plot(point_a, point_b, linewidth=linewidth,
                             zorder=0, alpha=0.1)

    def update_annotion(self, event, annot, agent_set):
        annot.xy = (event.xdata, event.ydata)
        info = self.agent_summary(agent_set)
        annot.set_text(info)

    def click(self, event, annot):
        vis = annot.get_visible()
        if event.inaxes == self.ax:
            for mpl_object, agent_set in self.mpl_collections.items():
                cont, ind = mpl_object.contains(event)
                if cont:
                    self.update_annotion(event, annot, agent_set)
                    annot.set_visible(True)
                    self.fig.canvas.draw_idle()
                    break
                else:
                    if vis:
                        annot.set_visible(False)
                        self.fig.canvas.draw_idle()

    def interactive_annotate(self):
        annot = self.ax.annotate('', xy=self.annot_property['xy'],
                                 xytext=self.annot_property['xytext'],
                                 textcoords=self.annot_property['textcoords'],
                                 bbox=self.annot_property['bbox'],
                                 verticalalignment='top')
        annot.set_visible(False)
        self.fig.canvas.mpl_connect('button_press_event',
                                    lambda event:
                                    self.click(event, annot))

    def plot_archetypes(self):
        mpl_collections = {}
        for i, arch in enumerate(self.archs):
            if arch == 'Reactor':
                reactors = self.plot_reactors()
                mpl_collections = {**mpl_collections, **reactors}
            else:
                non_reactors = self.plot_nonreactors(arch, i)
                mpl_collections = {**mpl_collections, **non_reactors}
        return mpl_collections

    def agent_summary(self, agent_set):
        # reactor marker for power output
        # transaction for transactions (commodity and avg amount)
        agent_set = sorted(list(agent_set))
        summary = ', '.join(agent_set)
        for i, agent in enumerate(agent_set):
            summary += "\n"
            agent = str(agent)
            name = self.agent_info[agent][0]
            spec = self.agent_info[agent][1]
            coordinates = (self.agent_info[agent]
                           [3], self.agent_info[agent][2])
            summary += "[" + spec + "]\n"
            summary += "Name: " + name + "\n"
            if spec == 'Reactor':
                capacity = str(self.reactor_markers[coordinates] /
                               self.capacity_to_markersize)
                summary += "Plant Capacity: " + capacity + " [MWe]\n"
        return summary


def main(sqlite_file):
    """ Calls all the functions above to produce the map output. Saves the
    resulting map as an html file.

    Parameters
    ----------
    sqlite_file: str
            path to cyclus output sqlite file

    Returns
    -------
    """
    cycvis = Cycvis(sqlite_file)
    cycvis.plot_transactions()
    cycvis.resize_legend()
    cycvis.interactive_annotate()
    plt.show()


if __name__ == "__main__":
    print('Running Directly')
    print('Usage: python cycvis.py [cyclus_output.sqlite]')
    main(sys.argv[1])
