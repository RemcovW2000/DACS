"""
Decorators for dependency tracking and lazy evaluation.
"""
from typing import Callable, Optional, Set, Any
from .base import TrackedInstance, _dependency_graph


class AttributeDescriptor:
    """Descriptor for tracked attributes with lazy evaluation."""

    def __init__(self, func: Callable, dependencies: Optional[Set[str]] = None):
        self.func = func
        self.name = func.__name__
        self.dependencies = dependencies or set()
        self._currently_computing = set()  # Track recursive computation

    def __set_name__(self, owner, name):
        self.name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self

        # Ensure instance is tracked
        if not isinstance(instance, TrackedInstance):
            raise TypeError(f"Attribute {self.name} can only be used on TrackedInstance subclasses")

        attr_id = instance._get_attr_id(self.name)

        # Check if we have a cached value
        if self.name in instance._cached_values:
            return instance._cached_values[self.name]

        # Prevent infinite recursion
        if attr_id in self._currently_computing:
            raise RecursionError(f"Circular dependency detected for attribute {self.name}")

        try:
            self._currently_computing.add(attr_id)

            # Track dependencies during computation
            old_tracking_context = getattr(_current_computation, 'attr_id', None)
            _current_computation.attr_id = attr_id

            # Compute the value
            value = self.func(instance)

            # Cache the result
            instance._cached_values[self.name] = value
            instance._attribute_versions[self.name] = instance._attribute_versions.get(self.name, 0) + 1

            return value

        finally:
            self._currently_computing.discard(attr_id)
            if 'old_tracking_context' in locals():
                _current_computation.attr_id = old_tracking_context
            else:
                _current_computation.attr_id = None

    def __set__(self, instance, value):
        """Setting an attribute manually invalidates cache and dependents."""
        if not isinstance(instance, TrackedInstance):
            raise TypeError(f"Attribute {self.name} can only be used on TrackedInstance subclasses")

        # Clear cached value
        if self.name in instance._cached_values:
            del instance._cached_values[self.name]

        # Store the value directly
        instance._tracked_attributes[self.name] = value

        # Invalidate dependents
        instance._invalidate_dependents(self.name)

    def __delete__(self, instance):
        """Delete cached value and stored value."""
        if self.name in instance._cached_values:
            del instance._cached_values[self.name]
        if self.name in instance._tracked_attributes:
            del instance._tracked_attributes[self.name]
        instance._invalidate_dependents(self.name)


class ComputationContext:
    """Context to track what attribute is currently being computed."""
    def __init__(self):
        self.attr_id = None

_current_computation = ComputationContext()


def Attribute(dependencies: Optional[Set[str]] = None):
    """
    Decorator for creating tracked attributes with lazy evaluation.

    Args:
        dependencies: Optional set of attribute names this attribute depends on.
                     If None, dependencies will be auto-detected during computation.

    Example:
        class MyClass(TrackedInstance):
            def __init__(self, x, y):
                super().__init__()
                self._x = x
                self._y = y

            @tracked_property
            def x(self):
                return self._x

            @x.setter
            def x(self, value):
                self._x = value

            @Attribute(dependencies={'x', 'y'})
            def sum(self):
                return self.x + self.y

            @Attribute()  # Auto-detect dependencies
            def doubled_sum(self):
                return 2 * self.sum  # Will auto-detect dependency on 'sum'
    """
    def decorator(func: Callable) -> AttributeDescriptor:
        return AttributeDescriptor(func, dependencies)
    return decorator


def track_access(instance: TrackedInstance, attr_name: str):
    """
    Track that the current computation accesses the given attribute.
    This should be called when an attribute is accessed during computation.
    """
    if not hasattr(_current_computation, 'attr_id') or _current_computation.attr_id is None:
        return

    current_attr_id = _current_computation.attr_id
    accessed_attr_id = instance._get_attr_id(attr_name)

    # Don't track self-dependency
    if current_attr_id != accessed_attr_id:
        _dependency_graph.add_dependency(current_attr_id, accessed_attr_id)


class TrackedProperty:
    """Property descriptor that automatically tracks access for dependency detection."""

    def __init__(self, fget=None, fset=None, fdel=None, doc=None):
        self.fget = fget
        self.fset = fset
        self.fdel = fdel
        self.__doc__ = doc

    def __set_name__(self, owner, name):
        self.name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self
        if self.fget is None:
            raise AttributeError(f"unreadable attribute {self.name}")

        # Track access for dependency detection
        if isinstance(instance, TrackedInstance):
            track_access(instance, self.name)

        return self.fget(instance)

    def __set__(self, instance, value):
        if self.fset is None:
            raise AttributeError(f"can't set attribute {self.name}")

        self.fset(instance, value)

        # Invalidate dependents if this is a tracked instance
        if isinstance(instance, TrackedInstance):
            instance._invalidate_dependents(self.name)

    def __delete__(self, instance):
        if self.fdel is None:
            raise AttributeError(f"can't delete attribute {self.name}")
        self.fdel(instance)

        # Invalidate dependents if this is a tracked instance
        if isinstance(instance, TrackedInstance):
            instance._invalidate_dependents(self.name)

    def getter(self, fget):
        return type(self)(fget, self.fset, self.fdel, self.__doc__)

    def setter(self, fset):
        return type(self)(self.fget, fset, self.fdel, self.__doc__)

    def deleter(self, fdel):
        return type(self)(self.fget, self.fset, fdel, self.__doc__)


def tracked_property(func):
    """
    Decorator to create a tracked property that automatically detects dependencies.
    Similar to @property but tracks access for dependency detection.
    """
    return TrackedProperty(func)


def tracked_attribute(initial_value: Any = None):
    """
    Simpler decorator that automatically creates a tracked property with getter and setter.

    Usage:
        class MyClass(TrackedInstance):
            width = tracked_attribute(5.0)  # Creates property with getter/setter
            height = tracked_attribute(3.0)

            @Attribute()
            def area(self):
                return self.width * self.height
    """
    class TrackedAttributeDescriptor:
        def __init__(self, default_value):
            self.default_value = default_value
            self.name = None
            self.private_name = None

        def __set_name__(self, owner, name):
            self.name = name
            self.private_name = f"_{name}"

        def __get__(self, instance, owner):
            if instance is None:
                return self

            # Track access for dependency detection
            if isinstance(instance, TrackedInstance):
                track_access(instance, self.name)

            # Return the value, using default if not set
            return getattr(instance, self.private_name, self.default_value)

        def __set__(self, instance, value):
            # Set the private attribute
            setattr(instance, self.private_name, value)

            # Invalidate dependents if this is a tracked instance
            if isinstance(instance, TrackedInstance):
                instance._invalidate_dependents(self.name)

    return TrackedAttributeDescriptor(initial_value)
