def discardinglogger(name, message):
    pass

def printinglogger(name, message):
    if name is None:
        print(message)
    else:
        print("%s: %s" % (name, message))

def getlogger(logger):
    if logger == "discard":
        return discardinglogger
    elif logger == "print" or logger is None:
        return printinglogger
    else:
        return logger
