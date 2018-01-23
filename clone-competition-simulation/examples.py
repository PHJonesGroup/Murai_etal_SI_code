from competition_and_mutation import Competition, MoranStyleComp, normal_fitness_dist, uniform_fitness_dist
from colourscales import get_colourscale_with_random_mutation_colour
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


def example1():
    # Run a single simulation of algorithm 1
    # If all parameters kept the same here, will be the first run from the multiple runs example
    np.random.seed(0)  # Change the seed or comment this line to see a different result

    WT_init = 9900  # Initial population of wildtype cells
    A_init = 100  # Initial population of cells with mutation A
    A_fitness = 1.05  # Initial fitness of cells with mutation A. Wildtype has fitness 1.
    mutation_rate = 0.3  # Chance of a single mutation per generation.
    fitness_dist = normal_fitness_dist(var=0.1, mean=1)  # The distribution for fitness of new mutations
    # fitness_dist = uniform_fitness_dist(0.8, 1.2)  # Uncomment to use a uniform distribution instead
    num_generations = 150  # Change to vary the length of the simulation
    combine_mutations = True  # Change to False to have new mutation fitness replace old instead of combining

    # Add a file path here to save the plot created. Use None to have plot pop up in window instead
    plot_file = 'example1_plot.pdf'

    c = Competition(WT_init=WT_init, A_init=A_init, A_fitness=A_fitness, mutation_rate=mutation_rate,
                    fitness_distribution=fitness_dist, num_generations=num_generations,
                    combine_mutations=combine_mutations, plot_file=plot_file,
                    figsize=(10, 8))
    c.run_sim()
    max_fitness = c.clones_array[:, c.growth_idx].max()
    cs = get_colourscale_with_random_mutation_colour(max_fitness)
    c.colourscales = cs
    c.plot()


def example2():
    # Run a single simulation of algorithm 2
    # This can be significantly slower than example1 above
    np.random.seed(2)  # Change the seed or comment this line to see a different result

    WT_init = 9900  # Initial population of wildtype cells
    A_init = 100  # Initial population of cells with mutation A
    A_fitness = 1.1  # Initial fitness of cells with mutation A. Wildtype has fitness 1.
    mutation_rate = 0.00025  # Chance of adding a mutation.
    fitness_dist = normal_fitness_dist(var=0.1, mean=1)  # The distribution for fitness of new mutations
    # fitness_dist = uniform_fitness_dist(0.8, 1.2)  # Uncomment to use a uniform distribution instead
    num_births = 1000000  # Change to vary the length of the simulation
    sampling = 500  # Change to vary how often the samples for the plot are taken. Small=more memory, Large=courser plot
    combine_mutations = True  # Change to False to have new mutation fitness replace old instead of combining

    plot_file = 'example2_plot.pdf'  # Add a file path here to save a pdf of the plot created

    c = MoranStyleComp(WT_init=WT_init, A_init=A_init, A_fitness=A_fitness, mutation_rate=mutation_rate,
                       fitness_distribution=fitness_dist, num_births=num_births, sampling=sampling,
                       combine_mutations=combine_mutations, plot_file=plot_file,
                       figsize=(10, 8))
    c.run_sim()
    max_fitness = c.clones_array[:, c.growth_idx].max()
    cs = get_colourscale_with_random_mutation_colour(max_fitness)
    c.colourscales = cs
    c.plot()


def multiple_runs_example():
    # Plotting multiple results on the same plot, the proportion of descendants of the labelled mutant population
    # This will output the plots to pdf files

    np.random.seed(0)  # Change the seed or comment this line to see a different result
    num_repeats = 12
    WT_init = 9900  # Initial population of wildtype cells
    A_init = 100  # Initial population of cells with mutation A
    A_fitness = 1.05  # Initial fitness of cells with mutation A. Wildtype has fitness 1.
    mutation_rate = 0.3  # Chance of a single mutation per generation.
    fitness_distribution = normal_fitness_dist(var=0.1, mean=1)  # The distribution for fitness of new mutations
    num_generations = 150  # Change to vary the length of the simulation
    combine_mutations = True  # Change to False to have new mutation fitness replace old instead of combining

    x = range(num_generations + 1)

    fig1, ax1 = plt.subplots()

    for n in range(num_repeats):
        c = Competition(WT_init=WT_init, A_init=A_init, A_fitness=A_fitness, mutation_rate=mutation_rate,
                        fitness_distribution=fitness_distribution, num_generations=num_generations,
                        combine_mutations=combine_mutations, figsize=(10, 8),
                        plot_file='single_plot{0}.pdf'.format(n))
        c.run_sim()
        max_fitness = c.clones_array[:, c.growth_idx].max()
        cs = get_colourscale_with_random_mutation_colour(max_fitness)
        c.colourscales = cs
        c.plot()
        A_proportion = c.proportional_populations[np.where(c.clones_array[:, c.type_idx] == 1)].sum(axis=0)
        if n == 0:
            ax1.plot(x, A_proportion, label=n, c='g', linewidth=3, zorder=2)
        else:
            ax1.plot(x, A_proportion, label=n, c=cm.Greys((n+2)/(num_repeats+5)), zorder=1)

    ax1.set_ylim([0, 1])
    ax1.set_xlim([0, num_generations])
    ax1.set_yticklabels([0, 20, 40, 60, 80, 100])
    ax1.set_xticks([])
    legend = ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig1.savefig('multiplot.pdf', bbox_extra_artists=(legend,), bbox_inches='tight')


if __name__ == '__main__':
    example1()
    # example2()
    # multiple_runs_example()

