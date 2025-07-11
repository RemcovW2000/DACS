import csv

# 1. Define your element groups, each with its name, direction, and element list
groups = [
    {
        "name": "TE top",
        "direction": "LE to TE",
        "elements": [i for i in range(18851, 22752, 100)]
    },
    {
        "name": "TE bottom",
        "direction": "LE to TE",
        "elements": [i for i in range(22850, 26751, 100)]
    },
    {
        "name": "Mid bottom",
        "direction": "TE to LE",
        "elements": [i for i in reversed(range(51, 3952, 100))]
    },
    {
        "name": "Mid top",
        "direction": "TE to LE",
        "elements": [i for i in reversed(range(8850, 12751, 100))]
    },
    {
        "name": "LE",
        "direction": "top to bottom",
        "elements": [i for i in range(18751, 12850, -100)]
    },
    {
        "name": "Front spar",
        "direction": "top to bottom",
        "elements": [i for i in range(8751, 6450, -100)]
    },
    {
        "name": "Rear spar",
        "direction": "top to bottom",
        "elements": [i for i in range(6350, 4049, -100)]
    },
]
