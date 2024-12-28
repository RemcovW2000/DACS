from abc import ABC, abstractmethod
from typing import Literal

failure_mode_options = Literal[
    'fiber_failure',
    'inter_fiber_failure',
    'first_ply_failure',
    'buckling',
    'wrinkling',
    'child'
]

name_options = Literal[
    'lamina',
    'laminate',
    'sandwich',
    'member',
    'airfoil',
    'wing',
]

failure_state_options = Literal[True, False]

class StructuralEntity(ABC):
    def __init__(self, name: name_options):
        # Failure indicator by default has 'child' for the potential child object(s) of the class, how can I add more
        # failure indicators for different failure modes of the object? this would be different for each class,
        # for example the member class would have a buckling failure mode, which I would like to add to this dictionary

        # furthermore, how would I overwrite  the failure indicator?
        self.failure_indicators = {
            'child': None,
        }
        self.name = name

    @property
    def child_objects(self):
        return []
    @abstractmethod
    def failure_analysis(self):
        '''
        Method that does a failure analysis of the parent class, and calls this method for all child classes.

        Is an augmentor class -> influences the state of the class

        Sets the correct failure indicator for the right failure mode(s)

        :return: None
        '''
        pass

    def set_failure_indicator(self, failure_mode: failure_mode_options, failure_indicator: float):
        # failure indicator is initialised as an empty dictionary, but then failure indicators are added to it
        self.failure_indicators[failure_mode] = failure_indicator

    def get_hierarchy(self):
        '''
        Returns the lower hierarchy of child objects
        :return:
        '''
        return