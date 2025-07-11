"""
Utility functions for dependency tracking system.
"""
from typing import Dict, Set, Any
from .base import _dependency_graph, TrackedInstance


def get_dependency_graph() -> Dict[str, Set[str]]:
    """Get a copy of the current dependency graph."""
    return {k: v.copy() for k, v in _dependency_graph.dependencies.items()}


def clear_dependency_graph():
    """Clear the entire dependency graph. Use with caution."""
    _dependency_graph.dependencies.clear()
    _dependency_graph.dependents.clear()


def print_dependencies(instance: TrackedInstance, attr_name: str = None):
    """Print dependency information for an instance or specific attribute."""
    if attr_name:
        attr_id = instance._get_attr_id(attr_name)
        deps = _dependency_graph.get_dependencies(attr_id)
        dependents = _dependency_graph.get_dependents(attr_id)

        print(f"Attribute: {attr_name}")
        print(f"  Depends on: {[aid.split('.')[-1] for aid in deps]}")
        print(f"  Dependents: {[aid.split('.')[-1] for aid in dependents]}")
    else:
        print(f"All dependencies for instance {instance.__class__.__name__}:")
        instance_prefix = f"{instance._instance_id}."

        for attr_id, dependents in _dependency_graph.dependencies.items():
            if attr_id.startswith(instance_prefix):
                attr_name = attr_id[len(instance_prefix):]
                filtered_dependents = [dep for dep in dependents if dep.startswith(instance_prefix)]
                if filtered_dependents:
                    dependent_names = [dep[len(instance_prefix):] for dep in filtered_dependents]
                    print(f"  {attr_name} -> {dependent_names}")


def invalidate_cache(instance: TrackedInstance, attr_name: str):
    """Manually invalidate the cache for a specific attribute."""
    if attr_name in instance._cached_values:
        del instance._cached_values[attr_name]
    instance._invalidate_dependents(attr_name)


def get_cached_attributes(instance: TrackedInstance) -> Dict[str, Any]:
    """Get all currently cached attribute values."""
    return instance._cached_values.copy()


def is_cached(instance: TrackedInstance, attr_name: str) -> bool:
    """Check if an attribute value is currently cached."""
    return attr_name in instance._cached_values
