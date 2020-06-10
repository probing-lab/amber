
from random import random


def run_program():
    i = 0
    x = 10
    while x > 0:
        i = i + 1
        if random() < 2/3:
            x = x + i
        else:
            x = x - i
        if x > 1000000000000:
            return False

    return True


count_term = 0
count_non_term = 0
for n in range(100):
    term = run_program()
    if term:
        count_term += 1
    else:
        count_non_term += 1
    print(f'{n} / {100}')


print(f'# Termination: {count_term}')
print(f'# Non-Termination: {count_non_term}')