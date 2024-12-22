from Toolbox.Laminate import LaminateBuilder

Laminate = LaminateBuilder([0,0,0,45], True, True, 1)

Laminate.Loads = [10, 0, 0, 0 ,0 ,0]

print(Laminate.FailureAnalysis())