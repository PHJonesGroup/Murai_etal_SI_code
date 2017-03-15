from competition_and_mutation import Competition, MoranStyleComp, normal_fitness_dist, uniform_fitness_dist
import matplotlib.pyplot as plt
import numpy as np


def example1():
    # Run a single simulation of algorithm 1
    np.random.seed(34)  # Change the seed or comment this line to see a different result

    WT_init = 9900  # Initial population of wildtype cells
    A_init = 100  # Initial population of cells with mutation A
    A_fitness = 1.05  # Initial fitness of cells with mutation A. Wildtype has fitness 1.
    mutation_rate = 0.2  # Chance of a single mutation per generation.
    fitness_dist = normal_fitness_dist(var=0.1, mean=1)  # The distribution for fitness of new mutations
    # fitness_dist = uniform_fitness_dist(0.8, 1.2)  # Uncomment to use a uniform distribution instead
    num_generations = 150  # Change to vary the length of the simulation
    combine_mutations = False  # Change to True to combine effects of multiple mutations

    plot_file = None  # Add a file path here to save a pdf of the plot created

    c = Competition(WT_init=WT_init, A_init=A_init, A_fitness=A_fitness, mutation_rate=mutation_rate,
                    fitness_distribution=fitness_dist, num_generations=num_generations,
                    combine_mutations=combine_mutations, plot_file=plot_file,
                    figsize=(10, 8))
    c.run_sim()
    c.plot()


def example2():
    # Run a single simulation of algorithm 2
    np.random.seed(0)  # Change the seed or comment this line to see a different result

    WT_init = 9900  # Initial population of wildtype cells
    A_init = 100  # Initial population of cells with mutation A
    A_fitness = 1.1  # Initial fitness of cells with mutation A. Wildtype has fitness 1.
    mutation_rate = 0.0001  # Chance of a single mutation per generation.
    fitness_dist = normal_fitness_dist(var=0.1, mean=1)  # The distribution for fitness of new mutations
    # fitness_dist = uniform_fitness_dist(0.8, 1.2)  # Uncomment to use a uniform distribution instead
    num_births = 1000000  # Change to vary the length of the simulation
    sampling = 100  # Change to vary how often the samples for the plot are taken. Small=more memory, Large=courser plot
    combine_mutations = False  # Change to True to combine effects of multiple mutations

    plot_file = None  # Add a file path here to save a pdf of the plot created

    c = MoranStyleComp(WT_init=WT_init, A_init=A_init, A_fitness=A_fitness, mutation_rate=mutation_rate,
                       fitness_distribution=fitness_dist, num_births=num_births, sampling=sampling,
                       combine_mutations=combine_mutations, plot_file=plot_file,
                       figsize=(10, 8))
    c.run_sim()
    c.plot()


def multiple_runs_example():
    # Plotting multiple results on the same plot, similar to Supplementary Figure 7a
    # This will output the plots to pdf files

    np.random.seed(1)  # Change the seed or comment this line to see a different result
    num_repeats = 10
    WT_init=9900
    A_init=100
    A_fitness=1.05
    mutation_rate=0.2
    fitness_distribution=normal_fitness_dist(0.1)
    num_generations=200
    combine_mutations=False

    x = xrange(num_generations + 1)

    fig1, ax1 = plt.subplots()

    for n in xrange(num_repeats):
        c = Competition(WT_init=WT_init, A_init=A_init, A_fitness=A_fitness, mutation_rate=mutation_rate,
                        fitness_distribution=fitness_distribution, num_generations=num_generations,
                        combine_mutations=combine_mutations, figsize=(10, 8),
                        plot_file='single_plot'+str(n))
        c.run_sim()
        c.plot()
        A_proportion = c.proportional_populations[np.where(c.clones_array[:, c.type_idx] == 1)].sum(axis=0)
        ax1.plot(x, A_proportion, label=n)

    ax1.set_ylim([0, 1])
    ax1.set_xlim([0, num_generations])
    legend = ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # fig1.show()
    fig1.savefig('multiplot.pdf', bbox_extra_artists=(legend,), bbox_inches='tight')


if __name__ == '__main__':
    example1()
    # example2()
    # multiple_runs_example()


