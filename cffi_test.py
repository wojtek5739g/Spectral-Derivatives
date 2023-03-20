import cffi_example

if __name__ == "__main__":
    x, y = 6, 2.3

    answer = cffi_example.lib.cmult(x, y)
    print(f' In Python: int: {x} float {y:.1f} return val {answer:.1f}')