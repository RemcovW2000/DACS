"""
Simple examples showing different ways to create tracked properties.
"""

from dtrack import TrackedInstance, Attribute, tracked_property, tracked_attribute, TrackedAttributes


# Method 1: The verbose way (what you saw before)
class BeamVerbose(TrackedInstance):
    def __init__(self, width, height):
        super().__init__()
        self._width = width
        self._height = height

    @tracked_property
    def width(self):
        return self._width

    @width.setter
    def width(self, value):
        self._width = value

    @tracked_property
    def height(self):
        return self._height

    @height.setter
    def height(self, value):
        self._height = value

    @Attribute()
    def area(self):
        print("Calculating area...")
        return self.width * self.height


# Method 2: Using tracked_attribute (much simpler!)
class BeamSimple(TrackedInstance):
    def __init__(self, width=5.0, height=3.0):
        super().__init__()
        # Set the initial values
        self.width = width
        self.height = height

    # These automatically create getters and setters with dependency tracking
    width = tracked_attribute(5.0)
    height = tracked_attribute(3.0)

    @Attribute()
    def area(self):
        return self.width * self.height

    @Attribute()
    def perimeter(self):
        return 2 * (self.width + self.height)


# Method 3: Using TrackedAttributes (group them together)
class BeamGrouped(TrackedInstance):
    def __init__(self, width=5.0, height=3.0, density=7850):
        super().__init__()
        # Set initial values
        self.props.width = width
        self.props.height = height
        self.props.density = density

    # Define all tracked properties in one place
    props = TrackedAttributes(
        width=5.0,
        height=3.0,
        density=7850
    )

    @Attribute()
    def area(self):
        return self.props.width * self.props.height

    @Attribute()
    def volume_per_length(self):
        return self.area  # area per unit length

    @Attribute()
    def mass_per_length(self):
        return self.volume_per_length * self.props.density


def demo_simple_usage():
    print("=== Simple Usage Demo ===")

    # All three approaches work the same way:
    beam1 = BeamVerbose(10, 5)
    beam2 = BeamSimple(10, 5)
    beam3 = BeamGrouped(10, 5, 8000)

    print("Initial areas:")
    print(f"Verbose beam: {beam1.area}")
    print(f"Simple beam: {beam2.area}")
    print(f"Grouped beam: {beam3.area}")

    print("\nChanging width to 20:")
    beam1.width = 20
    beam2.width = 20
    beam3.props.width = 20

    print(f"Verbose beam: {beam1.area}")
    print(f"Simple beam: {beam2.area}")
    print(f"Grouped beam: {beam3.area}")