"""
Example usage of the dependency tracking system.

This example demonstrates how to use the @Attribute decorator and TrackedInstance
to create a class with automatic dependency tracking and lazy evaluation.
"""

from dtrack import TrackedInstance, Attribute, tracked_property, print_dependencies


class Rectangle(TrackedInstance):
    """Example class demonstrating dependency tracking."""

    def __init__(self, width: float, height: float):
        super().__init__()
        self._width = width
        self._height = height

    @tracked_property
    def width(self):
        """Width property with automatic dependency tracking."""
        return self._width

    @width.setter
    def width(self, value):
        self._width = value

    @tracked_property
    def height(self):
        """Height property with automatic dependency tracking."""
        return self._height

    @height.setter
    def height(self, value):
        self._height = value

    @Attribute()
    def area(self):
        """Area computed lazily with auto-detected dependencies."""
        print(f"Computing area: {self.width} * {self.height}")
        return self.width * self.height

    @Attribute()
    def perimeter(self):
        """Perimeter computed lazily with auto-detected dependencies."""
        print(f"Computing perimeter: 2 * ({self.width} + {self.height})")
        return 2 * (self.width + self.height)

    @Attribute()
    def diagonal(self):
        """Diagonal computed lazily with auto-detected dependencies."""
        print(f"Computing diagonal: sqrt({self.width}² + {self.height}²)")
        return (self.width ** 2 + self.height ** 2) ** 0.5

    @Attribute()
    def area_to_perimeter_ratio(self):
        """Ratio computed lazily, depends on other computed attributes."""
        print(f"Computing ratio: {self.area} / {self.perimeter}")
        return self.area / self.perimeter


class Circle(TrackedInstance):
    """Another example class to show the system works with different classes."""

    def __init__(self, radius: float):
        super().__init__()
        self._radius = radius

    @tracked_property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        self._radius = value

    @Attribute()
    def area(self):
        """Circle area: π * r²"""
        import math
        print(f"Computing circle area: π * {self.radius}²")
        return math.pi * self.radius ** 2

    @Attribute()
    def circumference(self):
        """Circle circumference: 2 * π * r"""
        import math
        print(f"Computing circumference: 2 * π * {self.radius}")
        return 2 * math.pi * self.radius

    @Attribute()
    def diameter(self):
        """Circle diameter: 2 * r"""
        print(f"Computing diameter: 2 * {self.radius}")
        return 2 * self.radius


def demo_basic_usage():
    """Demonstrate basic usage of the dependency tracking system."""
    print("=== Basic Usage Demo ===")

    # Create a rectangle
    rect = Rectangle(5, 3)

    print("1. First access to area (will compute):")
    print(f"Area: {rect.area}")

    print("\n2. Second access to area (cached):")
    print(f"Area: {rect.area}")

    print("\n3. Access perimeter (will compute):")
    print(f"Perimeter: {rect.perimeter}")

    print("\n4. Access ratio (will compute, uses cached area and perimeter):")
    print(f"Ratio: {rect.area_to_perimeter_ratio}")

    print("\n5. Change width (invalidates dependent attributes):")
    rect.width = 10

    print("6. Access area again (will recompute due to width change):")
    print(f"Area: {rect.area}")

    print("7. Access ratio again (will recompute due to area change):")
    print(f"Ratio: {rect.area_to_perimeter_ratio}")


def demo_dependency_visualization():
    """Demonstrate dependency visualization."""
    print("\n=== Dependency Visualization Demo ===")

    rect = Rectangle(4, 6)

    # Access all computed attributes to establish dependencies
    _ = rect.area
    _ = rect.perimeter
    _ = rect.diagonal
    _ = rect.area_to_perimeter_ratio

    print("\nDependency graph for rectangle:")
    print_dependencies(rect)


def demo_multiple_instances():
    """Demonstrate that different instances have separate dependency tracking."""
    print("\n=== Multiple Instances Demo ===")

    rect1 = Rectangle(2, 3)
    rect2 = Rectangle(4, 5)
    circle = Circle(3)

    print("Rectangle 1 area:", rect1.area)
    print("Rectangle 2 area:", rect2.area)
    print("Circle area:", circle.area)

    print("\nChanging rectangle 1 width:")
    rect1.width = 10

    print("Rectangle 1 area (recomputed):", rect1.area)
    print("Rectangle 2 area (unchanged, cached):", rect2.area)
    print("Circle area (unchanged, cached):", circle.area)


if __name__ == "__main__":
    demo_basic_usage()
    demo_dependency_visualization()
    demo_multiple_instances()
