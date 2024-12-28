def addone(lst):
    # Check if list is empty:
    if not lst:
        raise ValueError("List is empty")

    # make new list to make sure we dont change the original list
    newlst = list(lst)
    all_nines = True
    for i, val in enumerate(newlst[::-1]):
        if val < 9:
            newlst[-(i+1)] += 1
            all_nines = False
            break
        else:
            newlst[-(i+1)] = 0

    if all_nines:
        newlst.insert(0, 1)
    return newlst
lst = [1,2,3,4,5]
print(addone(lst))
print(lst)
lst = [1,0,3,0,9]
print(addone(lst))
lst = [9,9,9,9]
print(addone(lst))
lst = []
print(addone(lst))

