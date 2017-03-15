from __future__ import division
import bisect
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class GeneralComp(object):
    """Base class for the simulations. Contains general functions included in both methods."""

    def __init__(self, WT_init, A_init, A_fitness, mutation_rate, fitness_distribution, combine_mutations,
                 sim_length, figsize=(10, 10), plot_file=None):
        """

        :param WT_init: Integer. The number of wildtype cells at the start of the simulation
        :param A_init: Integer. The number of cells containing mutation A at the start of the simulation
        :param A_fitness: Float>0. The "fitness" of mutation A. Has a slightly different interpretation in the two
         algorithms
        :param mutation_rate: Float between 0 and 1. Controls the rate at which new mutations are introduced.
         Again, slightly different usage in the two algorithms.
        :param fitness_distribution: A function which returns a random positive float when called.
        :param combine_mutations: Boolean. If True, the effects of mutations will be multiplied when in the same cell.
         If False, the fitness value of the last mutation added to a cell will replace the previous value
        :param sim_length: Integer. The length of the simulation.
        :param figsize: Tuple, size of the plot to output
        :param plot_file: None or String. Path to a file to save the plot as a pdf - .pdf will be added to the filename
         If None, will show the plot instead of saving it.
        """

        # Clones array  - stores information about each clone
        # Columns: Clone id  | Type (WT or A) | Fitness | Generation Born | Parent Population
        self.clones_array = np.array([
            [0, 0, 1, 0, -1],  # WT has fitness 1. If no other clones, homeostatic.
            [1, 1, A_fitness, 0, -1]
        ])
        # Include indices here for later use
        self.id_idx = 0
        self.type_idx = 1
        self.growth_idx = 2
        self.generation_born_idx = 3
        self.parent_idx = 4

        self.sim_length = sim_length

        self.population_array = np.zeros((2, self.sim_length + 1))  # Will store the population counts
        self.population_array[0, 0] = WT_init  # Initial WT pop
        self.population_array[1, 0] = A_init  # Initial mutation A pop
        self.initial_pop = WT_init + A_init
        self.A_init_proportion = A_init / self.initial_pop

        self.mutation_rate = mutation_rate  # Exact meaning depends on simulation method.
        self.fitness_distribution = fitness_distribution  # Function that returns a fitness for a new mutant

        # Whether the effects of multiple mutations are combined or if the last mutation replaces any previous
        self.combine_mutations = combine_mutations

        # Details for plotting
        self.figsize = figsize
        self.plot_order = []  # Used to generate the plot of clones
        self.descendant_counts = {}
        self.plot_file = plot_file

    def get_new_fitness(self, old_fitness):
        """

        :param old_fitness: Float. The previous fitness value in the cell which has acquired the new mutation.
         If self.combine_mutations is True, the effects of mutations will be multiplied.
        :return: Float
        """
        random_draw = self.fitness_distribution()
        if self.combine_mutations:
            new_fitness = old_fitness * random_draw
        else:
            new_fitness = random_draw
        return new_fitness

    # The functions for plotting the results
    def plot(self):
        self.split_populations_for_plot()  # Breaks up the populations so subclones appear from their parent clone
        fig, ax = plt.subplots(figsize=self.figsize)

        self.cumulative_array = np.cumsum(self.split_pops_for_plotting, axis=0)
        self.get_colours()  # Get the colours used for the plots

        self.make_stackplot(ax)  # Add the clone populations to the plots

        # Add the clone births to the plot
        x = []
        y = []
        c = []
        for clone in self.clones_array:
            gen = clone[self.generation_born_idx]
            if gen > 0:
                pops = np.where(self.plot_order == clone[self.id_idx])[0][0]
                x.append(gen)
                y.append(self.split_pops_for_plotting[:pops][:, int(gen)].sum())
                c.append('r')

        ax.scatter(x, y, c=c, marker='x')
        plt.ylim([0, 1])
        plt.xlim([0, self.sim_length])

        if self.plot_file:
            plt.savefig('{0}.pdf'.format(self.plot_file))
        else:
            plt.show()

    def get_colour(self, rate, clone_type):
        if clone_type == 0:
            return cm.YlOrBr(rate / 2)  # Colours for wildtype subclones
        else:
            return cm.Greens(rate)  # Colours for mutation A subclones

    def get_colours(self):
        # Generate the colours for the clones plot. Colour depends on type (wildtype/A) and relative fitness.
        rates = self.clones_array[:, self.growth_idx]
        min_ = rates.min() - 0.1
        max_ = rates.max()
        self.colours = {}
        for clone in self.clones_array:
            scaled_rate = (clone[self.growth_idx] - min_) / (max_ - min_)
            self.colours[clone[self.id_idx]] = self.get_colour(scaled_rate, clone[self.type_idx])

        return [self.colours[o] for o in self.plot_order]

    def order_results(self):
        # Order the arrays to put all mutation A, then all WT
        wt_clone_idx = self.clones_array[:, self.type_idx] == 0
        self.clones_array = np.concatenate((self.clones_array[~wt_clone_idx], self.clones_array[wt_clone_idx]),
                                           axis=0)
        self.population_array = np.concatenate(
            (self.population_array[~wt_clone_idx], self.population_array[wt_clone_idx]),
            axis=0)

    def get_children(self, idx):
        # Return the ids of subclones of the given clone
        return self.clones_array[self.clones_array[:, self.parent_idx] == idx][:, self.id_idx]

    def get_descendants(self, idx, order):
        # Find the subclones of the given clone. Runs iteratively until found all descendants
        # Adds clone ids to order list.
        # order will be used to make the stackplot so the subclones appear from their parent clone
        order.append(idx)
        children = self.get_children(idx)
        np.random.shuffle(children)
        self.descendant_counts[idx] = len(children)
        for ch in children:
            if ch != idx:
                self.get_descendants(ch, order)
                order.append(idx)

    def split_populations_for_plot(self):
        # Breaks up the populations so subclones appear from their parent clone
        # Mutation A clones first so they appear at the bottom of the plot
        order_A = []
        order_WT = []

        self.get_descendants(1, order_A)
        self.get_descendants(0, order_WT)

        A_split_array = np.concatenate([self.proportional_populations[self.clones_array[:, self.id_idx] == o] / \
                                          (self.descendant_counts[o] + 1) for o in order_A])
        WT_split_array = np.concatenate([self.proportional_populations[self.clones_array[:, self.id_idx] == o] / \
                                         (self.descendant_counts[o] + 1) for o in order_WT])
        self.split_pops_for_plotting = np.concatenate((A_split_array, WT_split_array), axis=0)

        self.plot_order = order_A + order_WT

    def make_stackplot(self, ax):
        # Make the stackplot using fill between
        for i in range(len(self.plot_order) - 1, -1, -1):
            colour = self.colours[self.plot_order[i]]
            array = self.cumulative_array[i]
            if i > 0:
                next_array = self.cumulative_array[i - 1]
            else:
                next_array = 0
            ax.fill_between(range(self.sim_length + 1), array, 0, where=array > next_array, facecolor=colour,
                            interpolate=True, linewidth=0)


class Competition(GeneralComp):
    """Runs a simulation of the clonal growth, mutation and competition. Algorithm 1"""

    def __init__(self, WT_init, A_init, A_fitness, mutation_rate, fitness_distribution,
                 combine_mutations, num_generations,
                 figsize=(10, 10), plot_file=None):
        """

        :param WT_init: Integer. The number of wildtype cells at the start of the simulation
        :param A_init: Integer. The number of cells containing mutation A at the start of the simulation
        :param A_fitness: Float>0. The fitness value for the cells containing mutation A at the start of the simulation.
         Fitness determines how much a clone will grow from one generation to the next. Wildtype cells have fitness 1.
        :param mutation_rate: Float between 0 and 1. The chance per generation that a single mutation will be added.
        :param fitness_distribution: A function which returns a random positive float when called.
        :param combine_mutations: Boolean. If True, the effects of mutations will be multiplied when in the same cell.
         If False, the fitness value of the last mutation added to a cell will replace the previous value
        :param num_generations: Integer. How many generations the simulation runs for.
        :param figsize: Tuple, size of the plot to output
        :param plot_file: None or String. Path to a file to save the plot as a pdf - .pdf will be added to the filename
         If None, will show the plot instead of saving it.
        """
        super(Competition, self).__init__(WT_init, A_init, A_fitness, mutation_rate, fitness_distribution,
                                          combine_mutations, sim_length=num_generations,
                                          figsize=figsize, plot_file=plot_file)
        self.generation = 0

    def run_sim(self):
        for i in range(self.sim_length):
            self.generation += 1

            # Calculate the next generation population sizes multiplying the previous population and the growth rates
            self.population_array[:, self.generation] = self.population_array[:,
                                                        self.generation - 1] * self.clones_array[:, self.growth_idx]

            # Normalise the population to match the original total population
            self.population_array[:, self.generation] = self.population_array[:, self.generation] / \
                                                        self.population_array[:,
                                                        self.generation].sum() * self.initial_pop

            # Remove populations once they drop below 1. Prevents any negative populations after a mutation
            self.population_array[:, self.generation][self.population_array[:, self.generation] < 1] = 0

            # Draw random number to decide if a mutation is introduced in this generation
            if np.random.random() < self.mutation_rate:
                self.introduce_mutation()

        # A steps to prepare for plotting
        self.order_results()
        self.proportional_populations = self.population_array / self.population_array.sum(axis=0)

    def introduce_mutation(self):
        """Select a fitness for the new mutation and the cell in which the mutation occurs"""
        # Randomly select a cell (equivalent to randomly selected a clone in proportion to the clone size)
        population_selector = np.random.random()
        cumsum = np.cumsum(self.population_array[:, self.generation], axis=0)
        idx = find_ge(cumsum, population_selector * cumsum[-1])
        selected_clone = self.clones_array[idx]
        selected_pop = self.population_array[idx]
        new_type = selected_clone[self.type_idx]  # Is the new clone labelled or not
        old_growth_rate = selected_clone[self.growth_idx]

        # Get a fitness value for the new clone.
        new_growth_rate = self.get_new_fitness(old_growth_rate)

        selected_pop[self.generation] -= 1  # Remove the mutated cell from the clone it was previously part of.
        new_clone = np.array([[len(self.clones_array), new_type, new_growth_rate,
                               self.generation, idx]])
        new_population = np.zeros((1, self.sim_length + 1))
        new_population[0, self.generation] = 1  # Start with a single cell population

        # Add the new clone to the arrays
        self.population_array = np.concatenate((self.population_array, new_population), axis=0)
        self.clones_array = np.concatenate((self.clones_array, new_clone), axis=0)


class MoranStyleComp(GeneralComp):
    """Runs a simulation of the clonal growth, mutation and competition. Algorithm 2"""

    def __init__(self, WT_init, A_init, A_fitness, mutation_rate, fitness_distribution,
                 combine_mutations, num_births, sampling,
                 figsize=(10, 10), plot_file=None):
        """

        :param WT_init: Integer. The number of wildtype cells at the start of the simulation
        :param A_init: Integer. The number of cells containing mutation A at the start of the simulation
        :param A_fitness: Float>0. The fitness value for the cells containing mutation A at the start of the simulation.
         Fitness determines how likely a cell is to be picked to replicate. Wildtype cells have fitness 1.
        :param mutation_rate: Float between 0 and 1.
         The chance that a mutation will be added before another birth-death.
         There can be multiple mutations between one birth-death and the next.
        :param fitness_distribution: A function which returns a random positive float when called.
        :param combine_mutations: Boolean. If True, the effects of mutations will be multiplied when in the same cell.
         If False, the fitness value of the last mutation added to a cell will replace the previous value
        :param num_births: Integer. The number of cell birth-deaths in the simulation.
        :param sampling: Integer. The number of steps simulated between each sample added to the plot.
         There are many more steps involved in this algorithm than in the other, as only one birth-death or mutation
         occurs per step. Retaining the populations after every step and plotting all of them would be too memory
         intensive.
        :param figsize: Tuple, size of the plot to output
        :param plot_file: None or String. Path to a file to save the plot as a pdf - .pdf will be added to the filename
         If None, will show the plot instead of saving it.
        """
        sim_length = int(num_births / sampling)
        super(MoranStyleComp, self).__init__(WT_init, A_init, A_fitness, mutation_rate, fitness_distribution,
                                             combine_mutations, sim_length, figsize=figsize, plot_file=plot_file)

        self.num_births = num_births  # Length of the simulation
        self.sampling = sampling  # How many births per sampling for the plots
        self.plot_idx = 0  # Keeping track of x-coordinate of the plot

    def run_sim(self):
        # Each step, pick one cell to be replicated and one to be replaced.
        # Replication chances depend on fitness
        # Replacement chances do not depend on fitness
        # E.g. for a population A with fitness r and i individuals out of N total, the probability of being selected
        # for reproduction is ri/(sum(r'*i')), where r' and i' are fitness and population of all clones
        # The probability of replacement is simply the proportion i/N

        # We assume that mutation is an externally driven event, so is independent of replication, replacement and
        # fitness.
        # This may be different from other similar processes where it is assumed there is a small chance of a mutation
        # at each cell division.
        # Between each birth-death, there is a small chance of a mutation occuring,
        # in which case a cell is selected at random to mutate and start a new clone.
        # The population array is then updated and we continue

        # The x-axis/time axis should not be taken very seriously.
        # It is merely a list of the order of events, and does not say anything about the speed of the results.

        # Much slower process in terms of iterations
        # Only single change per iteration
        # Instead of storing every single change, we only keep a record every n-steps, given by self.sampling
        current_population = self.population_array[:, 0]

        assert self.mutation_rate < 1, 'Mutation rate >= 1. While loop would never exit!'

        self.plot_idx = 0
        i = 0
        while i < self.num_births:
            if np.random.random() < self.mutation_rate:  # Random draw to determine if we introduce a new mutation
                current_population = self.introduce_mutation(current_population)
            else:
                i += 1
                # Select population to replicate cell
                # Select random number to select which population
                birth_selector = np.random.random()
                # make cumulative list of the fitnesses
                fitness_cumsum = np.cumsum(current_population * self.clones_array[:, self.growth_idx], axis=0)
                # Pick out the selected population
                birth_idx = find_ge(fitness_cumsum, birth_selector * fitness_cumsum[-1])

                # Select replaced population
                death_selector = np.random.random()
                cumsum = np.cumsum(current_population, axis=0)
                death_idx = find_ge(cumsum, death_selector * cumsum[-1])

                # Update population array
                current_population[birth_idx] += 1
                current_population[death_idx] -= 1

                if i % self.sampling == 0:  # Regularly take a sample for the plot
                    self.plot_idx += 1
                    self.population_array[:, self.plot_idx] = current_population

        # A steps to prepare for plotting
        self.order_results()
        self.proportional_populations = self.population_array / self.population_array.sum(axis=0)

    def introduce_mutation(self, current_population):
        """Select a fitness for the new mutation and the cell in which the mutation occurs"""
        # Randomly select a cell (equivalent to randomly selected a clone in proportion to the clone size)
        population_selector = np.random.random()
        cumsum = np.cumsum(current_population, axis=0)
        idx = find_ge(cumsum, population_selector * cumsum[-1])
        selected_clone = self.clones_array[idx]
        new_type = selected_clone[self.type_idx]  # Is the new clone labelled or not
        old_growth_rate = selected_clone[self.growth_idx]

        # Get a fitness value for the new clone.
        new_growth_rate = self.get_new_fitness(old_growth_rate)

        current_population[idx] -= 1  # Remove the mutated cell
        new_clone = np.array([[len(self.clones_array), new_type, new_growth_rate,
                               self.plot_idx, idx]])
        current_population = np.concatenate((current_population, [1]))  # Start with a single cell population

        # Add the new clone to the arrays
        self.population_array = np.concatenate((self.population_array,
                                                np.zeros((1, self.sim_length + 1))), axis=0)
        self.clones_array = np.concatenate((self.clones_array, new_clone), axis=0)
        return current_population


def find_ge(a, x):
    """Find leftmost item greater than or equal to x"""
    i = bisect.bisect_left(a, x)
    if i != len(a):
        return i
    raise ValueError


def normal_fitness_dist(var, mean=1):
    # Returns a function used to draw the new mutation rates from
    # Done this way so can be called without parameters
    # Allow us to replace with other distributions in the simulation if desired
    def fitness_func():
        g = np.random.normal(mean, var)
        if g < 0:
            print 'growth rate below zero! Redrawing a new rate'
            return fitness_func()
        return g

    return fitness_func


def uniform_fitness_dist(min_, max_):
    # Returns a function used to draw the new mutation rates from
    # Done this way so can be called without parameters
    # Allow us to replace with other distributions in the simulation if desired
    if min_ < 0:
        print 'Minimum of uniform distribution set below zero. Setting to zero instead.'
        min_ = 0

    def fitness_func():
        g = np.random.uniform(min_, max_)
        return g

    return fitness_func
