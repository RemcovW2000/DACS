class objects:
    def __init__(self, nr):
        object = None
        self.nr = nr

    def nothing(self):
        return

obj1 = objects(1)
obj2 = objects(2)
obj3 = objects(3)
def find_object_by_coordinate(objects_list, coordinate):
    for obj, end_coordinate in objects_list:
        if coordinate <= end_coordinate:
            return obj
    return None  # or raise an exception if you prefer

# Example usage
objects_list = [[obj1, 100], [obj2, 200], [obj3, 300]]

# Test the function
coordinate1 = 140
coordinate2 = 90

result1 = find_object_by_coordinate(objects_list, coordinate1).nr
result2 = find_object_by_coordinate(objects_list, coordinate2).nr

print(result1)  # Should print obj2
print(result2)  # Should print obj1