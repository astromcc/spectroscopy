def printhowdy(nickname):
    if type(nickname) == str:
        print("Howdy, " + nickname)
    else:
        print("Didn't catch that.")
    return

def printhello(nickname):
    if type(nickname) == str:
        print("Hello " + nickname)
    else:
        print("Didn't catch that.")
    return

print('simpleImport.py has run')