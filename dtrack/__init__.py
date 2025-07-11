"""
Dependency tracking system with lazy evaluation.

This package provides decorators and base classes for automatic dependency tracking
and lazy evaluation of computed attributes.
"""

from .base import TrackedInstance, _dependency_graph
from .decorators import Attribute, tracked_property, track_access, tracked_attribute
from .utils import (
    get_dependency_graph,
    clear_dependency_graph,
    print_dependencies,
    invalidate_cache,
    get_cached_attributes,
    is_cached
)

__all__ = [
    'TrackedInstance',
    'Attribute',
    'tracked_property',
    'tracked_attribute',
    'track_access',
    'get_dependency_graph',
    'clear_dependency_graph',
    'print_dependencies',
    'invalidate_cache',
    'get_cached_attributes',
    'is_cached'
]
