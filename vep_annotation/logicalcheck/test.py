

class Test:

    def __init__(self, name, age):
        self.name = name
        self.age = age 




if __name__ == '__main__':
    name = 'mike'
    age = 22

    print(Test(age, name).__dict__)