import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from Laminate import Laminate
from Lamina import Lamina

class MonteCarloSimulation:
    def __init__(self, num_simulations, mean_properties, std_dev_properties, laminas_configuration):
        """
        Initializes the Monte Carlo simulation with the number of simulations, mean properties,
        standard deviation properties, and the configuration of laminas.

        :param num_simulations: int, number of simulations to run
        :param mean_properties: dict, mean values of the laminate properties
        :param std_dev_properties: dict, standard deviations of the laminate properties
        :param laminas_configuration: list of tuples, each tuple contains the thickness, angle,
                                      and a dict with keys 'elastic' and 'failure' for properties
        """
        self.num_simulations = num_simulations
        self.mean_properties = mean_properties
        self.std_dev_properties = std_dev_properties
        self.laminas_configuration = laminas_configuration
        self.results = []

    def sample_properties(self):
        """
        Samples properties based on the specified mean and standard deviation.

        :return: dict, sampled properties
        """
        sampled_properties = {}
        for key in self.mean_properties:
            mean = self.mean_properties[key]
            std_dev = self.std_dev_properties[key]
            sampled_properties[key] = np.random.normal(mean, std_dev)
        return sampled_properties

    def run_simulation(self, load):
        """
        Runs a single simulation of the Monte Carlo analysis.

        :param load: numpy array, the load applied to the laminate
        """
        sampled_properties = self.sample_properties()
        laminas = []

        for thickness, angle, properties in self.laminas_configuration:
            # Ensure these property keys match those in the mean_properties and std_dev_properties dictionaries
            elastic_properties = [sampled_properties[prop] for prop in properties['elastic']]
            failure_properties = [sampled_properties[prop] for prop in properties['failure']]
            lamina = Lamina(thickness, angle, elastic_properties, failure_properties)
            laminas.append(lamina)

        laminate = Laminate(laminas) # Creating Laminate object
        laminate.Loads = load # Assigning load to the laminate

        # Performing failure analysis and storing the maximum failure factor
        failure_state, failed_laminas, maxfailfactor = laminate.FailureAnalysis()
        self.results.append(maxfailfactor)  # Storing the result of the simulation

    def run(self, load):
        """
        Runs the Monte Carlo simulations.

        :param load: numpy array, the load applied to the laminate
        """
        for i in tqdm(range(self.num_simulations)):
            # print('sim:', i)
            self.run_simulation(load)

    def analyze_results(self):
        """
        Analyzes and reports the results of the Monte Carlo simulations.
        """
        # Counting how many simulations resulted in at least one failure
        failures = sum(1 for result in self.results if result > 1)
        print(f"Failures: {failures}/{self.num_simulations}")
        # Calculating and printing the failure rate
        failure_rate = failures/self.num_simulations
        print('Failure rate: ', failure_rate*100, '%')
        return failures