from pyswip import Prolog


def queryProlog(query, consultList):
    prolog = Prolog()

    for file in consultList:
        prolog.consult(file)

    solutions = prolog.query(query)

    solutionlist = []
    for solution in solutions:
        solutionlist.append(solution)

    return solutionlist
