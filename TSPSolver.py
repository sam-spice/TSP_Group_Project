#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import copy



class TSPSolver:
    def __init__( self, gui_view ):
        self._scenario = None
        self.best_path_cost = None
        self.number_pruned = 0
        self.max_states = 0
        self.bssf_updates = 0

    def setupWithScenario( self, scenario ):
        self._scenario = scenario


    ''' <summary>
        This is the entry point for the default solver
        which just finds a valid random tour
        </summary>
        <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (
not counting initial BSSF estimate)</returns> '''
    def defaultRandomTour( self, start_time, time_allowance=60.0 ):

        results = {}


        start_time = time.time()

        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        while not foundTour:
            # create a random permutation
            perm = np.random.permutation( ncities )

            #for i in range( ncities ):
                #swap = i
                #while swap == i:
                    #swap = np.random.randint(ncities)
                #temp = perm[i]
                #perm[i] = perm[swap]
                #perm[swap] = temp

            route = []

            # Now build the route using the random permutation
            for i in range( ncities ):
                route.append( cities[ perm[i] ] )

            bssf = TSPSolution(route)
            #bssf_cost = bssf.cost()
            #count++;
            count += 1

            #if costOfBssf() < float('inf'):
            if bssf.costOfRoute() < np.inf:
                # Found a valid route
                foundTour = True
        #} while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
        #timer.Stop();

        results['cost'] = bssf.costOfRoute() #costOfBssf().ToString();                          // load results array
        results['time'] = time.time() - start_time
        results['count'] = count
        results['soln'] = bssf

       # return results;
        return results



    def greedy( self, start_time, time_allowance=60.0 ):
        pass

    def b_b_set_up(self):
        self.max_states = 0
        self.number_pruned = 0
        self.bssf_updates = 0
        cities = self._scenario.getCities()
        cities_array = [[0 for x in range(len(cities) + 1)] for y in range(len(cities) + 1)]
        for i in range(len(cities)):
            cities_array[i + 1][0] = cities[i]._name
            cities_array[0][i + 1] = cities[i]._name

        for i in range(len(cities)):
            for j in range(len(cities)):
                cities_array[i + 1][j + 1] = cities[i].costTo(cities[j])
        initial_state = BBstate()
        initial_state.cities_array = cities_array
        return initial_state

    def reduce(self, to_reduce):
        if len(to_reduce.cities_array) <= 3:
            return
        for row in to_reduce.cities_array[1:]:
            reduce_by = min(row[1:])
            to_reduce.cost += reduce_by
            for i in range(1,len(to_reduce.cities_array)):
                row[i] -= reduce_by

        for x in range(1, len(to_reduce.cities_array)):
            col = []
            for y in range(1, len(to_reduce.cities_array)):
                col.append(to_reduce.cities_array[y][x])
            reduce_by = min(col)
            to_reduce.cost += reduce_by
            for y in range(1, len(to_reduce.cities_array)):
                to_reduce.cities_array[y][x] -= reduce_by

    def branch_round(self, queue, state):
        for i in range(1, len(state.cities_array[0])):
            if state.current_city >= len(state.cities_array):
                continue
            cost = state.cities_array[state.current_city][i]
            if cost + state.cost > self.best_path_cost:
                self.number_pruned +=1
                continue

            new_row_char = state.cities_array[0][i]
            new_state = copy.deepcopy(state)
            new_state.cost += cost
            new_state.visited.append(new_row_char)
            old_row = new_state.cities_array[new_state.current_city][0]


            for j in range(len(new_state.cities_array)):
                new_state.cities_array[j].pop(i)

            new_state.cities_array.pop(new_state.current_city)
            new_state.current_city = math.inf
            for j in range(1, len(new_state.cities_array)):
                if new_state.cities_array[j][0] == new_row_char:
                    new_state.current_city = j

            for j in range(1, len(new_state.cities_array[0])):
                if new_state.cities_array[0][j] == old_row:
                    new_state.cities_array[new_state.current_city][j] = math.inf


            self.reduce(new_state)
            new_state.number_of_rounds += 1
            if new_state.cost <= self.best_path_cost:
                heapq.heappush(queue,new_state)
            else:
                self.number_pruned += 1

    def get_solution(self, visited):
        to_return = list()
        for entry in visited:
            for i in self._scenario.getCities():
                if entry == i._name:
                    to_return.append(i)
        return to_return


    def branchAndBound(self, start_time, time_allowance=60.0 ):
        bssf_results = self.defaultRandomTour(start_time, time_allowance)
        self.best_path_cost = bssf_results['cost']
        solution_count = 0
        bssf = bssf_results['soln']

        t_1 = time.time()
        initial_state = self.b_b_set_up()
        self.reduce(initial_state)
        pq = list()

        for i in range(1, len(initial_state.cities_array)):
            to_add = copy.deepcopy(initial_state)
            to_add.current_city = i
            to_add.visited.append(to_add.cities_array[i][0])
            pq.append(to_add)
        heapq.heapify(pq)

        while pq.__len__() > 0:
            popped = heapq.heappop(pq)
            if len(popped.cities_array) == 1:
                if popped.cost < self.best_path_cost:
                    self.bssf_updates += 1
                    cities = self.get_solution(popped.visited)
                    bssf = TSPSolution(cities[:-1])
                    self.best_path_cost = popped.cost
                solution_count += 1
            elif popped.cost > self.best_path_cost:
                self.number_pruned += 1
            else:
                self.branch_round(pq,popped)
            heapq.heapify(pq)
            self.max_states = max(self.max_states,len(pq))
            if time.time() - t_1 > time_allowance:
                break

        results = dict()
        results['cost'] = bssf.costOfRoute()
        t_2 = time.time()
        results['time'] = t_2 - t_1
        results['count'] = solution_count
        results['soln'] = bssf

        print('Max stored states: ' + str(self.max_states))
        print('pruned states: ' + str(self.number_pruned))
        print('bssf updates: ' + str(self.bssf_updates))
        return results

    def fancy(self, start_time, time_allowance=60.0):
        pass


class BBstate:
    def __init__(self):
        self.cities_array = []
        self.visited = []
        self.route = []
        self.cost = 0
        self.number_of_rounds = 0
        self.current_city = None
    def __lt__(self, other):
        if self.number_of_rounds == other.number_of_rounds:
            return self.cost < other.cost
        else:
            return self.number_of_rounds > other.number_of_rounds