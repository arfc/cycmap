from collections import defaultdict, OrderedDict
from mpl_toolkits.basemap import Basemap
import sqlite3 as sql
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import sys
import warnings
warnings.filterwarnings('ignore')


class Cycvis():
    # default settings for mpl plotting
    figsize = (11, 5.5)
    ax_main_bounds = [0.02, 0.05, 0.6, 0.8]
    ax_sub_bounds = [2, 2, 0, 0]
    ax_check_bounds = [0.02, 0.85, 0.6, 0.1]
    show_incommod = True
    show_outcommod = True

    def __init__(self, file_name):
        self.cur = self.get_cursor(file_name)
        self.agent_info = self.get_agent_info()
        self.archs = self.list_available_archetypes()
        self.init_yr, self.timestep, self.timestep_yr = self.get_sim_info()
        (self.fig, self.ax_main, self.ax_sub,
         self.ax_checkbox, self.agent_description) = self.setup_empty_plot()
        self.transactions = self.get_transactions()
        self.reactor_info = self.get_reactor_info()
        self.ax_main_plot_basemap()
        self.setup_checkbox()
        self.mpl_objects = self.ax_main_plot_agents()

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

    def list_available_archetypes(self):
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

    def get_sim_info(self):
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
        timestep_yr = init_year + (timestep / 12)
        return init_year, timestep, timestep_yr

    def get_arch_info(self, archetype):
        archetype_info = {}
        for k, v in self.agent_info.items():
            agentid = k
            spec = v[1]
            if spec == archetype:
                archetype_info[agentid] = v
        return archetype_info

    def calculate_basemap_bounds(self):
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
        bound_addition = 0.05
        latitudes = [info[2] for info in list(self.agent_info.values())]
        longitudes = [info[3] for info in list(self.agent_info.values())]
        llcrnrlat = min(latitudes)
        llcrnrlon = min(longitudes)
        urcrnrlat = max(latitudes)
        urcrnrlon = max(longitudes)
        extra_lon = (urcrnrlon - llcrnrlon) * bound_addition
        extra_lat = (urcrnrlat - llcrnrlat) * 3 * bound_addition
        llcrnrlat -= extra_lat
        llcrnrlon -= extra_lon
        urcrnrlat += extra_lat
        urcrnrlon += extra_lon
        bounds = [llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon]
        return bounds

    def get_reactor_info(self):
        query = ("SELECT DISTINCT timeseriespower.agentid, value,"
                 " longitude, latitude"
                 " FROM timeseriespower"
                 " INNER JOIN agentposition"
                 " ON timeseriespower.agentid = agentposition.agentid"
                 " WHERE value > 0")
        results = self.cur.execute(query)
        reactor_info = {}
        for row in results:
            reactor_info[row['agentid']] = [row['value'],
                                            (row['longitude'],
                                             row['latitude'])]
        return reactor_info

    def calculate_reactor_marker_size(self):
        """ Returns a dictionary of reactor coordinates and marker size using its
        power capacity

        Parameters
        ----------
        cur: sqlite cursor
                sqlite cursor

        Returns
        -------
        reactor_markers: dict
                dictionary with "
                key = (longitude, latitude) in float, and
                value = marker size in float for matplotlib scatter"
        """
        capacity_to_markersize = 0.5
        reactor_markers = defaultdict(float)
        for k, v in self.reactor_info.items():
            coordinates = v[1]
            capacity = v[0]
            reactor_markers[coordinates] += (capacity *
                                             capacity_to_markersize)
        return reactor_markers

    def get_transactions(self):
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
        query = ("SELECT senderid, receiverid, commodity, quantity, time"
                 " FROM TRANSACTIONS"
                 " INNER JOIN RESOURCES"
                 " ON TRANSACTIONS.resourceid = RESOURCES.resourceid"
                 " ORDER BY time")
        results = self.cur.execute(query)
        transaction_dict = defaultdict(list)
        for row in results:
            lifetime = self.agent_info[str(row['senderid'])][-1]
            if lifetime == -1:
                lifetime = (self.timestep[-1] -
                            self.agent_info[str(row['senderid'])][-2])
            transaction_dict[(row['senderid'],
                              row['receiverid'],
                              row['commodity']
                              )
                             ].append((row['time'], row['quantity']))
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
        arch_info = self.get_arch_info(arch)
        # agent[3] = longitude
        # agent[2] = latitude
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

    def calculate_transactions_line_width(self, arch):
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
        lines: dict
                dictionary of transactions and average quantity moved during
                lifetime with "
                key = tuple of sender coord, receiver coord, and commodity, and
                value = average quantity moved during lifetime in [MTHM]"
        """
        kg_to_tons = 1 / 1000
        lons, lats, labels = self.get_lons_lats_labels(arch, True)
        lines = defaultdict(float)
        for key in self.transactions.keys():
            commodity = key[2]
            sender_id = str(key[0])
            receiver_id = str(key[1])
            sender_lon = self.agent_info[sender_id][3]
            sender_lat = self.agent_info[sender_id][2]
            receiver_lon = self.agent_info[receiver_id][3]
            receiver_lat = self.agent_info[receiver_id][2]
            quantity = (np.sum([v[1] for v in self.transactions[key]]) *
                        kg_to_tons)
            sender_coord = (sender_lon, sender_lat)
            receiver_coord = (receiver_lon, receiver_lat)
            lines[(sender_coord, receiver_coord, commodity)] = quantity
        return lines

    def setup_empty_plot(self):
        fig = plt.figure(figsize=self.figsize)
        ax_main = fig.add_axes(self.ax_main_bounds)
        ax_sub = fig.add_axes(self.ax_sub_bounds)
        ax_checkbox = fig.add_axes(self.ax_check_bounds)
        annotation_bbox = dict(boxstyle="round", alpha=(0.0), fc="w")
        annotation = ax_main.annotate('',
                                      xy=(0, 0),
                                      visible=False,
                                      xytext=(1.02, 0.99),
                                      bbox=annotation_bbox,
                                      verticalalignment='top',
                                      textcoords='axes fraction')
        return fig, ax_main, ax_sub, ax_checkbox, annotation

    def setup_checkbox(self):
        # Set x, y ticker first
        self.ax_checkbox.get_xaxis().set_visible(False)
        self.ax_checkbox.get_yaxis().set_visible(False)
        # Set axis limit
        self.ax_checkbox.set_xlim(-0.05, 2)
        self.ax_checkbox.set_ylim(-0.5, 0.5)
        options = ['In Commod', 'Out Commod']
        x_text_offset = 0.18
        y_text_offset = -0.02
        for i, option in enumerate(options):
            self.ax_checkbox.scatter(i, 0)
            self.ax_checkbox.annotate(option,
                                      xy=(i + x_text_offset, y_text_offset),
                                      verticalalignment='center',
                                      horizontalalignment='center')

    def ax_main_plot_basemap(self):
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

        bounds = self.calculate_basemap_bounds()
        basemap = Basemap(ax=self.ax_main,
                          projection='cyl',
                          llcrnrlat=bounds[0],
                          llcrnrlon=bounds[1],
                          urcrnrlat=bounds[2],
                          urcrnrlon=bounds[3],
                          fix_aspect=False,
                          anchor='NW')
        basemap.drawstates(zorder=0)
        basemap.drawcountries(zorder=0)
        basemap.drawcoastlines(zorder=-15)
        basemap.fillcontinents(zorder=-5,
                               color='white',
                               lake_color='lightblue')
        basemap.drawmapboundary(zorder=-10,
                                fill_color='lightblue')

    def format_legend(self, ax):
        """ Resizes scatter plot legends to the same size

        Parameters
        ----------

        Returns
        -------
        """
        handles, labels = ax.get_legend_handles_labels()
        legend = OrderedDict(zip(labels, handles))
        legend = ax.legend(legend.values(),
                           legend.keys(),
                           columnspacing=0.1,
                           labelspacing=0.1,
                           fontsize='x-small',
                           loc='best')
        for handle in legend.legendHandles:
            handle._sizes = [30]
            handle._alpha = 1
            handle._linewidth = 3
        ax.figure.canvas.update()

    def ax_main_plot_reactors(self):
        """ Scatter plot of reactors with reactor capacity as marker size

        Parameters
        ----------

        Returns
        -------
        """
        mpl = {}
        lons, lats, labels = self.get_lons_lats_labels('Reactor', True)
        reactor_markers = self.calculate_reactor_marker_size()
        for i, (lon, lat, label) in enumerate(zip(lons, lats, labels)):
            agents = set(label.split(', '))
            mpl[self.ax_main.scatter(lon, lat,
                                     zorder=5,
                                     alpha=0.4,
                                     color='grey',
                                     label='Reactor',
                                     edgecolors='black',
                                     s=reactor_markers[(lon, lat)])] = agents
            self.ax_main.text(lon, lat, label,
                              fontsize=8,
                              verticalalignment='center',
                              horizontalalignment='center')
        return mpl

    def ax_main_plot_nonreactors(self, arch, i):
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
        num_archs = len(self.archs)
        colors = cm.rainbow(np.linspace(0, 1, num_archs))
        lons, lats, labels = self.get_lons_lats_labels(arch, True)
        for j, (lon, lat, label) in enumerate(zip(lons, lats, labels)):
            agents = set(label.split(', '))
            label = str(arch)
            mpl[self.ax_main.scatter(lon, lat,
                                     s=150,
                                     zorder=5,
                                     alpha=0.4,
                                     label=label,
                                     color=colors[i])] = agents
            self.ax_main.text(lon, lat, label,
                              fontsize=8,
                              verticalalignment='center',
                              horizontalalignment='center')
        return mpl

    def ax_main_plot_transactions(self):
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
        quantity_to_linewidth = 0.02
        commods = {k[2] for k in self.transactions.keys()}
        num_commods = len(commods)
        commod_cm = cm.Dark2(np.linspace(0, 1, num_commods))
        commod_cm = {commod: color for commod, color
                     in zip(commods, commod_cm)}
        for arch in self.archs:
            transactions = self.calculate_transactions_line_width(arch)
            for key, value in transactions.items():
                commod = key[2]
                point_a = [key[0][0], key[1][0]]
                point_b = [key[0][1], key[1][1]]
                linewidth = np.log(value * quantity_to_linewidth)
                self.ax_main.plot(point_a, point_b,
                                  zorder=0,
                                  alpha=0.1,
                                  label=commod,
                                  linewidth=linewidth,
                                  color=commod_cm[commod])

    def ax_sub_plot_agent_info(self, event, agent_set):
        self.agent_description.xy = (event.xdata, event.ydata)
        summary = self.ax_sub_get_agent_info(agent_set)
        self.agent_description.set_text(summary)
        self.agent_description.set_visible(True)

    def ax_sub_get_agent_info(self, agent_set):
        agent_set = sorted(list(agent_set))
        summary = ''
        for i, agent in enumerate(agent_set):
            name = self.agent_info[agent][0]
            spec = self.agent_info[agent][1]
            summary += "(" + agent + ") "
            summary += "[" + spec + "] "
            summary += name + "    "
            if spec == 'Reactor':
                capacity = str(self.reactor_info[int(agent)][0])
                summary += "(" + capacity + " [MWe])"
            summary += "\n"
        return summary

    def ax_main_click_event(self, event):
        visible = self.agent_description.get_visible()
        if event.inaxes == self.ax_main:
            for mpl_object, agent_set in self.mpl_objects.items():
                cont, ind = mpl_object.contains(event)
                if cont:
                    self.ax_sub_plot_agent_info(event, agent_set)
                    self.ax_sub_plot_transactions(agent_set)
                    self.fig.canvas.draw_idle()
                    break
                else:
                    if visible:
                        self.agent_description.set_visible(False)
                        self.ax_sub_update_bounds(self.ax_sub_bounds)
                        self.fig.canvas.draw_idle()

    def interactive_annotate(self):
        self.fig.canvas.mpl_connect('button_press_event',
                                    lambda event:
                                    self.ax_main_click_event(event))

    def ax_main_plot_agents(self):
        mpl_objects = {}
        for i, arch in enumerate(self.archs):
            if arch == 'Reactor':
                reactors = self.ax_main_plot_reactors()
                mpl_objects = {**mpl_objects, **reactors}
            else:
                non_reactors = self.ax_main_plot_nonreactors(arch, i)
                mpl_objects = {**mpl_objects, **non_reactors}
        return mpl_objects

    def calculate_cumulative_timeseries(self, in_list):
        """ returns a timeseries list from in_list data.

        Parameters
        ----------
        in_list: list
            list of data to be created into timeseries
            list[0] = time
            list[1] = value, quantity
        multiplyby: int
            integer to multiply the value in the list by for
            unit conversion from kilograms

        Returns
        -------
        timeseries of commodities in kg or tons
        """
        value = 0
        value_timeseries = []
        array = np.array(in_list)
        duration = 1 + int(self.timestep[-1] - self.timestep[0])
        for i in range(0, duration):
            value += sum(array[array[:, 0] == i][:, 1])
            value_timeseries.append(value)
        return value_timeseries

    def calculate_ax_sub_bounds(self):
        # https://stackoverflow.com/questions/
        # 29702424/how-to-get-matplotlib-figure-size
        self.ax_main.figure.canvas.draw()
        annot_box = self.agent_description.get_bbox_patch()
        annot_box_height = annot_box.get_height()
        invert_dpi_scale = self.fig.dpi_scale_trans.inverted()
        fig_bbox = self.fig.get_window_extent().transformed(invert_dpi_scale)
        fig_height = self.fig.dpi * fig_bbox.height
        annot_box_height_fraction = annot_box_height / fig_height
        margin = 0.05
        main_plot_width = self.ax_main_bounds[2] + self.ax_main_bounds[0]
        main_plot_height = self.ax_main_bounds[1] + self.ax_main_bounds[3]
        left = main_plot_width + margin
        bottom = 2 * margin
        width = 1 - (margin + left)
        height = main_plot_height - (annot_box_height_fraction +
                                     bottom +
                                     margin)
        bounds = [left, bottom, width, height]
        return bounds

    def ax_sub_get_plot_title(self, agent_set):
        title = []
        if self.show_incommod:
            title += ['In']
        if self.show_outcommod:
            title += ['Out']
        title = ' and '.join(title)
        title += ' Commodities of ' + ', '.join(agent_set)
        title += ' vs. Time'
        return title

    def ax_sub_update_bounds(self, new_bounds):
        self.ax_sub.remove()
        self.ax_sub = self.fig.add_axes(new_bounds)
        self.ax_sub.set_xlabel('Time [yr]', fontsize='small')
        self.ax_sub.set_ylabel('Mass [kg]', fontsize='small')
        self.ax_sub.tick_params(axis='both',
                                labelsize='x-small')
        self.ax_sub.ticklabel_format(axis='y',
                                     style='sci',
                                     useOffset=True,
                                     scilimits=(0, 0))
        self.ax_sub.yaxis.offsetText.set_fontsize('x-small')

    def ax_sub_get_label(self, commod, agentid, is_incommod):
        label = "(" + agentid + ") "
        if is_incommod:
            label += "In:  " + commod
        else:
            label += "Out: " + commod
        return label

    def ax_sub_plot_transactions(self, agent_set):
        sub_ax_coords = self.calculate_ax_sub_bounds()
        self.ax_sub_update_bounds(sub_ax_coords)
        for k, v in self.transactions.items():
            for agent in agent_set:
                senderid = str(k[0])
                receiverid = str(k[1])
                commod = k[2]
                if self.show_incommod and agent == receiverid:
                    is_in = True
                    label = self.ax_sub_get_label(commod, agent, is_in)
                    commod_timeseries = self.calculate_cumulative_timeseries(v)
                    self.ax_sub.plot(self.timestep_yr, commod_timeseries,
                                     label=label)
                if self.show_outcommod and agent == senderid:
                    is_in = False
                    label = self.ax_sub_get_label(commod, agent, is_in)
                    commod_timeseries = self.calculate_cumulative_timeseries(v)
                    self.ax_sub.plot(self.timestep_yr, commod_timeseries,
                                     label=label)
        self.ax_sub.set_title(self.ax_sub_get_plot_title(agent_set),
                              fontsize='small')
        self.format_legend(self.ax_sub)


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
    cycvis.ax_main_plot_transactions()
    cycvis.format_legend(cycvis.ax_main)
    cycvis.interactive_annotate()
    plt.show()


if __name__ == "__main__":
    print('Running Directly')
    print('Usage: python cycvis.py [cyclus_output.sqlite]')
    main(sys.argv[1])
