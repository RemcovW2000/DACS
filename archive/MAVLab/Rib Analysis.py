from Toolbox.laminate import laminate_builder

Laminate = laminate_builder([0,0,0,45], True, True, 1)

Laminate.Loads = [10, 0, 0, 0 ,0 ,0]

print(Laminate.failure_analysis())