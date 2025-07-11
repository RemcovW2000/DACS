"""
Base classes for dependency tracking system.
"""
from typing import Any, Set, Dict, Callable, Optional
import weakref
from functools import wraps


class DependencyGraph:
    """Global dependency graph to track relationships between attributes."""

    def __init__(self):
        self.dependencies: Dict[str, Set[str]] = {}  # attr_id -> set of dependent attr_ids
        self.dependents: Dict[str, Set[str]] = {}    # attr_id -> set of attributes it depends on

    def add_dependency(self, dependent: str, dependency: str):
        """Add a dependency relationship: dependent depends on dependency."""
        if dependency not in self.dependencies:
            self.dependencies[dependency] = set()
        if dependent not in self.dependents:
            self.dependents[dependent] = set()

        self.dependencies[dependency].add(dependent)
        self.dependents[dependent].add(dependency)

    def remove_dependency(self, dependent: str, dependency: str):
        """Remove a dependency relationship."""
        if dependency in self.dependencies:
            self.dependencies[dependency].discard(dependent)
        if dependent in self.dependents:
            self.dependents[dependent].discard(dependency)

    def get_dependents(self, attr_id: str) -> Set[str]:
        """Get all attributes that depend on the given attribute."""
        return self.dependencies.get(attr_id, set()).copy()

    def get_dependencies(self, attr_id: str) -> Set[str]:
        """Get all attributes that the given attribute depends on."""
        return self.dependents.get(attr_id, set()).copy()

    def clear_dependents(self, attr_id: str):
        """Clear all dependents of an attribute."""
        if attr_id in self.dependencies:
            for dependent in self.dependencies[attr_id].copy():
                self.remove_dependency(dependent, attr_id)


# Global dependency graph
_dependency_graph = DependencyGraph()


class TrackedInstance:
    """Base class for instances that support dependency tracking."""

    def __init__(self):
        self._tracked_attributes: Dict[str, Any] = {}
        self._attribute_versions: Dict[str, int] = {}
        self._cached_values: Dict[str, Any] = {}
        self._instance_id = id(self)

    def _get_attr_id(self, attr_name: str) -> str:
        """Get unique identifier for an attribute."""
        return f"{self._instance_id}.{attr_name}"

    def _invalidate_dependents(self, attr_name: str):
        """Invalidate all attributes that depend on the given attribute."""
        attr_id = self._get_attr_id(attr_name)
        dependents = _dependency_graph.get_dependents(attr_id)

        for dependent_id in dependents:
            # Extract instance_id and attr_name from dependent_id
            instance_id_str, dep_attr_name = dependent_id.rsplit('.', 1)
            instance_id = int(instance_id_str)

            # Find the instance (this is a simplified approach)
            # In a real implementation, you might want to use a registry
            if instance_id == self._instance_id:
                if dep_attr_name in self._cached_values:
                    del self._cached_values[dep_attr_name]
                if dep_attr_name in self._attribute_versions:
                    self._attribute_versions[dep_attr_name] += 1
