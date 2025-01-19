import math
import sys
import numpy as np

# Define the substituent structure


class SubStr:
    def __init__(self, g=0, m=0, sub_b=0, sub_m=0, o=' '):
        self.g = g          # Total number of carbons
        self.m = m          # Main chain length
        self.sub_b = sub_b  # Substituent balance
        self.sub_m = sub_m  # Substituent mass
        self.o = o          # Representation method (character)


# Global variable initialization
# Corresponds to C++'s ofstream fout("output.txt");
fout = open("output.txt", "w", encoding="utf-8")

Place = np.zeros([2, 201], dtype=int)


def power(expon, base=19):
    a = 1
    for _ in range(expon):
        a *= base
    return a


# Initialize substituent array s[1] to s[8]
sub_max = 16
s = [SubStr() for _ in range(20)]  # s[0] is unused, s[1] to s[9] are valid
s[0] = SubStr(g=0, m=0, sub_b=0, sub_m=0, o=' ')
s[1] = SubStr(g=1, m=1, sub_b=0, sub_m=0, o='C')
s[2] = SubStr(g=2, m=2, sub_b=0, sub_m=0, o='E')
s[3] = SubStr(g=3, m=2, sub_b=power(1), sub_m=1, o='I')
s[4] = SubStr(g=3, m=3, sub_b=0, sub_m=0, o='P')
s[5] = SubStr(g=4, m=2, sub_b=power(1) + power(2), sub_m=1, o='X')
s[6] = SubStr(g=4, m=3, sub_b=power(1), sub_m=1, o='y')
s[7] = SubStr(g=4, m=3, sub_b=power(3), sub_m=1, o='Y')
s[8] = SubStr(g=4, m=4, sub_b=0, sub_m=0, o='B')

s[9] = SubStr(g=5, m=3, sub_b=power(1) + power(2), sub_m=1, o='5')
s[10] = SubStr(g=5, m=3, sub_b=power(3) + power(4), sub_m=1, o='5')
s[11] = SubStr(g=5, m=3, sub_b=power(1) + power(3), sub_m=1, o='5')
s[12] = SubStr(g=5, m=3, sub_b=power(1) * 2, sub_m=2, o='5')
s[13] = SubStr(g=5, m=4, sub_b=power(1), sub_m=1, o='5')
s[14] = SubStr(g=5, m=4, sub_b=power(3), sub_m=1, o='5')
s[15] = SubStr(g=5, m=4, sub_b=power(5), sub_m=1, o='5')
s[16] = SubStr(g=5, m=5, sub_b=power(5), sub_m=1, o='5')


# Function to balance the structure
def balance(i, j, bal_, middle_odd, s_type, m_length):
    # center
    if j == middle_odd and m_length % 2 == 1:
        return bal_

    # left
    if j <= middle_odd:
        bal_ += power(2 * j - 1 + i) * s_type

    # right
    if j > middle_odd:
        pos = m_length - j + 1
        bal_ -= power(2 * pos - 1 + i) * s_type

    return bal_


def output(m_length, display_mod):
    if display_mod == 1:
        input("")
    elif display_mod != 2:
        return

    board = [[' ' for _ in range(80)] for _ in range(5)]

    for i in range(0, 2):
        for j in range(1, m_length + 1):
            if s[Place[i][j]].g > 0:
                if i == 0:
                    board[0][j * 3 - 2] = s[Place[i][j]].o
                    board[1][j * 3 - 2] = '|'
                else:
                    board[4][j * 3 - 2] = s[Place[i][j]].o
                    board[3][j * 3 - 2] = '|'
            board[2][j * 3 - 2] = 'C'
            board[2][j * 3 - 1] = '-'
            board[2][j * 3] = '-'

    for j in range(0, 5):
        line = ''.join(board[j][:m_length * 3 - 1])
        print(line)
        fout.write(line + '\n')
    print()
    fout.write('\n')


# Recursive function 'put' to check isomers
def put(last_s_type, last_position, left_carbon, bal, m_length, display_mod):
    counter = 0
    bal_ = bal

    # When left_carbon is 0, it means we've reached the end of the recursion and start checking the isomers
    if left_carbon == 0:
        if bal_ < 0:
            return 0

        # Check if the main chain should be swapped
        for i in range(2, m_length):
            # Left-side balance check
            bal_ = 0
            if i == s[Place[0][i]].m + 1 or i == s[Place[1][i]].m + 1:
                root = i
                for j in range(2, root):
                    bal_ += power(2 * (root - j) - 1) * Place[0][j]
                    bal_ += power(2 * (root - j)) * Place[1][j]

                if i == s[Place[0][i]].m + 1 and bal_ < s[Place[0][i]].sub_b:
                    return 0
                if i == s[Place[1][i]].m + 1 and bal_ < s[Place[1][i]].sub_b:
                    return 0

            # Right-side balance check
            bal_ = 0
            if m_length - i == s[Place[0][i]].m or m_length - i == s[Place[1][i]].m:
                root = i
                for j in range(root + 1, m_length):
                    bal_ += power(2 * (j - root) - 1) * Place[0][j]
                    bal_ += power(2 * (j - root)) * Place[1][j]

                if m_length - i == s[Place[0][i]].m and bal_ < s[Place[0][i]].sub_b:
                    return 0
                if m_length - i == s[Place[1][i]].m and bal_ < s[Place[1][i]].sub_b:
                    return 0

        output(m_length, display_mod)
        return 1

    else:
        middle_odd = (m_length + 1) // 2

        if not last_position:  # First substituent chain
            for s_type in range(sub_max, 0, -1):  # From sub_max to 1
                if s[s_type].g <= left_carbon:  # If there are remaining carbons
                    for j in range(2, m_length):
                        if j > s[s_type].m and m_length - j + 1 > s[s_type].m:

                            bal_new = balance(
                                0, j, bal_, middle_odd, s_type, m_length)
                            Place[0][j] = s_type
                            counter += put(s_type, (0, j), left_carbon -
                                           s[s_type].g, bal_new, m_length, display_mod)
                            Place[0][j] = 0

        else:        # For other substituent chains
            for s_type in range(last_s_type, 0, -1):
                if s[s_type].g <= left_carbon:
                    for i in range(0, 2):
                        for j in range(2, m_length):

                            # The first row must have the corresponding second row to be placed
                            if i == 1 and Place[0][j] == 0:
                                continue

                            # The current position must be empty
                            if not (Place[i][j] == 0):
                                continue

                            # Either the current substituent type must be different from the original, or the placement must be after it
                            if s_type == last_s_type and (i * m_length + j) < (last_position[0] * m_length + last_position[1]):
                                continue

                            # Cannot exceed the length of the main chain
                            if j <= s[s_type].m or m_length - j + 1 <= s[s_type].m:
                                continue

                            bal_new = balance(
                                i, j, bal_, middle_odd, s_type, m_length)

                            Place[i][j] = s_type
                            counter += put(s_type, (i, j), left_carbon -
                                           s[s_type].g, bal_new, m_length, display_mod)
                            Place[i][j] = 0

        return counter

# 主函数


def main():
    print("This program simulates the isomers of alkanes.")
    print("About substituents:")
    print(" Ethyl      Displayed as E")
    print(" n-Propyl   Displayed as P")
    print(" Isopropyl  Displayed as I")
    print(" n-Butyl    Displayed as X")
    print(" t-Butyl    Displayed as y")
    print(" Iso-butyl  Displayed as Y")
    print(" Butyl      Displayed as B")

    # "A: 1 Find isomers of a specified alkane, 2 Find isomers up to a specified number of carbons"
    # "B: 1 Display each isomer and wait, 2 Display all isomers, 3 Show the number of isomers for each main chain length, 4 Display only the number of isomers"
    mod_a = 1
    display_mod = 2

    def main_loop():
        counter = 0
        b = 0
        print(f"For C{tot_max}H{tot_max * 2 + 2}, there are isomers:")
        fout.write(f"For C{tot_max}H{tot_max * 2 + 2}, there are isomers:\n")
        for i in range(tot_max, 0, -1):
            m_length = i
            b = put(0, None, tot_max - m_length, 0, m_length, display_mod)
            counter += b
            if display_mod == 3:
                input("")
            print(f"When the main chain has {m_length} carbons, there are {b} isomers. So far, {counter} isomers have been found.")
            fout.write(f"When the main chain has {m_length} carbons, there are {b} isomers. So far, {counter} isomers have been found.\n\n")
        print(f"There are {counter} isomers in total.")
        fout.write(
            f"There are {counter} isomers in total.\n========================\n\n\n")

    while True:
        tot_max_input = input(
            "Please enter the total number of carbons for the alkane to be calculated:")
        try:
            tot_max = int(tot_max_input)
        except ValueError:
            print("Please enter a valid number.")
            continue

        if mod_a == 1:
            # Mode A=1: Find isomers of a specified alkane
            main_loop()

        elif mod_a == 2:
            # Mode A=2: Find isomers up to a specified number of carbons
            total_max_ = tot_max
            tot_max = 1
            while tot_max <= total_max_:
                main_loop()
                tot_max += 1

        continue_input = input(
            "Do you want to continue calculating? (y/n):").strip().lower()
        if continue_input != 'y':
            break

    fout.close()


if __name__ == "__main__":
    main()
