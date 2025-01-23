def square(num:int)->int:
    if num ==0:
        return 0
    return num**2 + square(num-1)

def main()->None:
    num = int(input("Enter the number \n"))
    print(f"The sum of squares if {square(num)}")

if __name__ == "__main__":
    main()
