def testReturn(inputbox1, inputbox2):
    text = ""
    text = text + inputbox2 + ": \n"
    for i in range(10):
        text = text + inputbox1 + '\n'
    print(text)
    return text