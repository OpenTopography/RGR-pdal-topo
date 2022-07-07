'''
Implements a priority queue using heapq. Python has a priority queue module built in, but it
is not stabley sorted (meaning that two items who are tied in priority are treated arbitrarily, as opposed to being
returned on a first in first out basis). This circumvents that by keeping a count on the items inserted and using
that count as a secondary priority
'''

import heapq

class priorityQueue:

    def __init__(self):
        '''
        Implementation of a priority queue that relies on heapqueue to create a stabley sorted priority queue. This
        ensures that items with equal priority will be 'popped' in the order that they were added to the queue.

        Each item in the heapqueue is a list of three items

        The first is the actual priorty, the second is a count of when this item was inserted, the third is the item.

        Heapqueue sorts based on these items successively, so when popping from the list the priority is considered,
        and then the count order.
        '''

        self.__pqueue = []
        self.__insertion_count = 0
        self.__nItems = 0

    def get(self):
        '''
        Remove an item and its priority from the queue and return both
        :return: (priority, item) - the highest priority item in the queue and its priority
        '''

        priority, count, item = heapq.heappop(self.__pqueue)
        self.__nItems -= 1
        return priority, item

    def put(self, priority, item):
        '''
        Add item to the priority queue with the specified priority

        :param priority: the priority with which this item should be added (low numbers = high priority)
        :param item: the item to add into the priority queue
        :return:
        '''
        self.__insertion_count += 1
        self.__nItems += 1
        entry = [priority, self.__insertion_count, item]
        heapq.heappush(self.__pqueue, entry)

    def isEmpty(self):
        '''
        Checks if the priority queue has an items remaining in it.
        :return: isEmpty: boolean, True if all items have been removed from the priority queue.
        '''
        return self.__nItems == 0