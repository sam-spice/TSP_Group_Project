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
        self.state_created = 0

    def setupWithScenario( self, scenario ):
        self._scenario = scenario


    ''' <summary>
        This is the entry point for the default solver
        which just finds a valid random tour
        </summary>
        <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (
not counting initial BSSF estimate)</returns> '''

    def defaultRandomTour(self, start_time, time_allowance=60.0):

        results = {}

        start_time = time.time()

        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        while not foundTour:

            if time.time() > (start_time + time_allowance):
                results['cost'] = np.inf
                results['soln'] = None
                results['count'] = 0
                results['time'] = time.time() - start_time
                return

            # create a random permutation
            perm = np.random.permutation(ncities)
            to_visit = np.ones(ncities)
            # perm = np.arange(ncities)

            # print('\n\nNEW PERMUTATION: {}'.format(perm))
            failed = False
            for i in range(ncities):
                # print('In city {}'.format(i))
                src = perm[i]
                dest = perm[(i + 1) % ncities]
                if self._scenario._edge_exists[src, dest]:
                    to_visit[dest] = False
                    continue
                elif i + 1 == ncities:
                    # print('CAN\'T GET BACK TO START!')
                    failed = True  # can't get back to start so try a new permutation
                    break

                else:  # try to swap the destination with a reachable one
                    # print('EDGE {}-->{} DOESN\'T EXIST'.format(src,dest))
                    reachable = np.nonzero(self._scenario._edge_exists[src, :] * to_visit)[0]
                    # print(reachable)
                    if np.sum(reachable > 0) == 0:
                        # print('DEAD END')
                        failed = True  # can't get back to start so try a new permutation
                        break
                    swapind = reachable[np.random.randint(np.sum(reachable > 0))]
                    # print('BEFORE: {}'.format(perm))
                    perm_loc_of_swapind = np.nonzero(perm == swapind)
                    perm[(i + 1) % ncities] = perm[perm_loc_of_swapind]
                    perm[perm_loc_of_swapind] = dest
                    to_visit[swapind] = False
                    # to_visit[swapind] = False
                    # print('AFTER: {}'.format(perm))
                    # print('REACHABLE: {}, picked {}'.format(reachable,swapind))

            if failed:
                # print('Trying a new permutation')
                continue  # try a new permutation

            route = []

            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])

            bssf = TSPSolution(route)
            # bssf_cost = bssf.cost()
            # count++;
            count += 1

            # if costOfBssf() < float('inf'):
            if bssf.costOfRoute() < np.inf:
                # Found a valid route
                foundTour = True
        # } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
        # timer.Stop();

        results['cost'] = bssf.costOfRoute()  # costOfBssf().ToString();                          // load results array
        results['time'] = time.time() - start_time
        results['count'] = count
        results['soln'] = bssf

        # return results;
        return results

    def greedy_setup(self):
        cities_list = copy.deepcopy(self._scenario.getCities())
        return cities_list


    def greedy_round(self, cities, start_city_idx):
        start_city = cities[start_city_idx]
        cities.remove(start_city)
        bssf_route = list()
        bssf_route.append(start_city)
        curr_city = start_city

        while(len(cities) > 0):
            min_dist = math.inf
            min_city = None
            for city in cities:
                if curr_city.costTo(city) <= min_dist:
                    min_dist = curr_city.costTo(city)
                    min_city = city
            bssf_route.append(min_city)
            cities.remove(min_city)
            curr_city = min_city
        #bssf_route.append(start_city)
        bssf = TSPSolution(bssf_route)
        return bssf

    def greedy( self, start_time, time_allowance=60.0 ):
        t_1 = time.time()
        bssf_cost = math.inf
        bssf = None
        count = 0
        cities_list = self.greedy_setup()
        for i in range(len(cities_list)):
            bssf_result = self.greedy_round(copy.deepcopy(cities_list),i)
            cost_of_route = bssf_result.costOfRoute()
            if bssf_result.costOfRoute() < bssf_cost:
                bssf_cost = bssf_result.costOfRoute()
                bssf = bssf_result
                count += 1
            if time.time() - t_1 > time_allowance:
                break
        t_2 = time.time()
        results = dict()
        results['cost'] = bssf.costOfRoute()  # costOfBssf().ToString();                          // load results array
        results['time'] = t_2 - t_1
        results['count'] = count
        results['soln'] = bssf
        return results




    # function called to set up the initial state array
    # time: O(n^2)
    # space: O(n^2)
    def b_b_set_up(self):
        # initializing all of the counter variables to 0
        self.max_states = 0
        self.number_pruned = 0
        self.bssf_updates = 0
        self.state_created = 0

        # get the list of cities and create the n x n array
        # time: O(n^2)
        # space O(n^2)
        cities = self._scenario.getCities()
        cities_array = [[0 for x in range(len(cities) + 1)] for y in range(len(cities) + 1)]

        # label the cities along the border
        # time: O(n)
        for i in range(len(cities)):
            cities_array[i + 1][0] = cities[i]._name
            cities_array[0][i + 1] = cities[i]._name

        # fill in the cost array
        # time: O(n^2)
        for i in range(len(cities)):
            for j in range(len(cities)):
                cities_array[i + 1][j + 1] = cities[i].costTo(cities[j])

        # instantiate the state object
        initial_state = BBstate()
        initial_state.cities_array = cities_array
        return initial_state

    # reduces the given state matrix and adjusts the cost accordingly
    # time: O(n^2)
    # space: O(1)
    def reduce(self, to_reduce):
        #if the array is small enough not worth reducing
        if len(to_reduce.cities_array) <= 3:
            return
        # reduce by rows
        # time: (n^2)
        # space: O(1)
        for row in to_reduce.cities_array[1:]:
            reduce_by = min(row[1:])
            to_reduce.cost += reduce_by
            for i in range(1,len(to_reduce.cities_array)):
                row[i] -= reduce_by

        # reduce by columns
        # time: O(n^2)
        # space O(1)
        for x in range(1, len(to_reduce.cities_array)):
            col = []
            for y in range(1, len(to_reduce.cities_array)):
                col.append(to_reduce.cities_array[y][x])
            reduce_by = min(col)
            to_reduce.cost += reduce_by
            for y in range(1, len(to_reduce.cities_array)):
                to_reduce.cities_array[y][x] -= reduce_by

    # expands the appropriate state and either adds the expansions to the queue or prunes them
    # time: O(n^3)
    # space: O(n^3)
    def branch_round(self, queue, state):

        # loops across the columns creating a new state for each column
        # representing travel from current city to other cities
        # time: O(n^3)
        # space: O(n^3)
        for i in range(1, len(state.cities_array[0])):
            # checks to see if the current city is out of bounds, should never happen but never hurts to check
            if state.current_city >= len(state.cities_array):
                continue
            self.state_created += 1
            # checks to see if the cost of the transfer + the current cost exceeds the BSSF, if it does prune it yo!
            cost = state.cities_array[state.current_city][i]
            if cost + state.cost > self.best_path_cost:
                self.number_pruned +=1
                continue

            # get the char representing the city
            new_row_char = state.cities_array[0][i]
            # create new copy of te state to be updated with the correct information!
            # space: O(n^2) because it needs to copy the entire array down which is at worst n^2
            # time: O(n^2) because it needs to copy the entire array
            new_state = copy.deepcopy(state)

            # updates the cost, and the row visited
            new_state.cost += cost
            new_state.visited.append(new_row_char)

            # stores the character of the old row for deletion purposes later
            old_row = new_state.cities_array[new_state.current_city][0]

            # gets rids of the entire column that is being visited next
            # time: O(n)
            # space: O(1)
            for j in range(len(new_state.cities_array)):
                new_state.cities_array[j].pop(i)

            # remove row jus used
            new_state.cities_array.pop(new_state.current_city)
            # set current row to infinite to make sure that it doesn't try to run for a city that no longer exists
            new_state.current_city = math.inf

            # find the row corresponding to the the next city to be visited and store it in the state object
            # time: O(n)
            for j in range(1, len(new_state.cities_array)):
                if new_state.cities_array[j][0] == new_row_char:
                    new_state.current_city = j

            # runs through to set the path from the next city back to the current city to be unusable
            # time: O(n)
            for j in range(1, len(new_state.cities_array[0])):
                if new_state.cities_array[0][j] == old_row:
                    new_state.cities_array[new_state.current_city][j] = math.inf

            # reduce newly created state
            # time: O(n^2)
            self.reduce(new_state)

            # increase number of rounds so that the priority queue functions correctly
            new_state.number_of_rounds += 1
            # checks if the new state should be pruned or added to the queue
            # time: O(logn)
            if new_state.cost <= self.best_path_cost:
                heapq.heappush(queue,new_state)
            else:
                self.number_pruned += 1

    # converts the characters being stored to actual city objects once the route is complete
    # time: O(n^2)
    # space: O(n)
    def get_solution(self, visited):
        to_return = list()
        # loops across city list and finds the cities in the order they are visited
        # time: O(n^2)
        # space: O(n)
        for entry in visited:
            for i in self._scenario.getCities():
                if entry == i._name:
                    to_return.append(i)
        return to_return

    # function that runs the branch and bound algorithm
    # time: O(n! * n^3)
    # space: O(n! * n^3)
    def branchAndBound(self, start_time, time_allowance=60.0 ):
        # run the default random tour algorithm to get baseline bssf for pruning purposes
        bssf_results = self.defaultRandomTour(start_time, time_allowance)
        self.best_path_cost = bssf_results['cost']
        solution_count = 0
        bssf = bssf_results['soln']

        # get initial time and then set up the initial state
        t_1 = time.time()
        initial_state = self.b_b_set_up()
        self.reduce(initial_state)
        # initialie the list to be used as a priority queue
        pq = list()

        # expands the initial state, one state for each possible starting city
        # time: O(n^3)
        # space: O(n^3)
        for i in range(1, len(initial_state.cities_array)):
            # time: O(n^2)
            # space: O(n^2)
            to_add = copy.deepcopy(initial_state)
            to_add.current_city = i
            to_add.visited.append(to_add.cities_array[i][0])
            pq.append(to_add)
            self.state_created += 1
        # time: O(n) I am gonna assume they do it in linear time because the python peeps a smart!
        heapq.heapify(pq)

        # run the reduction loop until the
        # time: O(n! * n^3) there are n! possible combos, as such it will take that many permutations * worst case time
        # space: O(n! * n^2) max number of combos * size of states
        while pq.__len__() > 0:
            # pop next item off the queue
            # time: O(logn)
            popped = heapq.heappop(pq)
            # if the city state has completed a tour
            # time: O(n^3) worst case for if statements
            # space: O(n^3)
            if len(popped.cities_array) == 1:
                if popped.cost < self.best_path_cost:
                    self.bssf_updates += 1
                    # time: O(n^2)
                    # space: O(n)
                    cities = self.get_solution(popped.visited)
                    bssf = TSPSolution(cities[:-1])
                    self.best_path_cost = popped.cost
                else:
                    self.number_pruned += 1
                solution_count += 1
            # is the state now needing to be pruned?
            elif popped.cost > self.best_path_cost:
                self.number_pruned += 1
            # incomplete state does not need to be pruned
            else:
                # time: O(n^3)
                # space: O(n^3)
                self.branch_round(pq,popped)

            #heapq.heapify(pq)
            # check if the number of states now in the queue is larger than the recorded maximum
            self.max_states = max(self.max_states,len(pq))
            # checks if time has elapsed if it has update the number of pruned to include all of the items in the list
            if time.time() - t_1 > time_allowance:
                self.number_pruned += len(pq)
                break


        # return the required info in a dictionary
        results = dict()
        results['cost'] = bssf.costOfRoute()
        t_2 = time.time()
        results['time'] = t_2 - t_1
        results['count'] = solution_count
        results['soln'] = bssf

        # prints more info related to the running of the algorithm to the terminal
        print('Max stored states: ' + str(self.max_states))
        print('pruned states: ' + str(self.number_pruned))
        print('bssf updates: ' + str(self.bssf_updates))
        print('total states created ' +  str(self.state_created)+'\n\n')

        return results

    def fancy(self, start_time, time_allowance=60.0):
        pass


class BBstate:
    def __init__(self):
        self.cities_array = []
        self.visited = []
        self.cost = 0
        self.number_of_rounds = 0
        self.current_city = None
    def __lt__(self, other):
        if self.number_of_rounds == other.number_of_rounds:
            return self.cost < other.cost
        else:
            return self.number_of_rounds > other.number_of_rounds