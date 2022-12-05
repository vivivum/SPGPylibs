class cexcept(Exception):
    def __init__(self, message, *args, **kwargs):
        super().__init__(message)

        for key, value in kwargs.items():
            setattr(self, key, value)

        for key, value in self.__dict__.items():
            print(key, value).format("%s => %s") 

def test_ce():
    message = "Exception Triggered! Something went wrong."
    raise cexcept(message)

def triggerException(num):
    if (num == 0):
        message = "Exception Triggered! Something went wrong."
        raise cexcept(message)
    else:
        print(num)


try:
    triggerException(0)
    print("Code has successfully been executed.")
except cexcept:
    print("Error: Number should not be 0.")