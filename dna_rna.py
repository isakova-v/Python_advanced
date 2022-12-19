from abc import ABC, abstractmethod
import numpy as np
from typing import Final

RNA_compare: Final[dict] = {'A': 'T', 'U': 'A', 'G': 'C', 'C': 'G'}
DNA_compare: Final[dict] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def is_dna(inp_str):
    inp_str1 = ""
    inp_str2 = ""
    try:
        inp_str1 = str(inp_str[0])
        inp_str2 = str(inp_str[1])
        inp_str = np.array(inp_str)
    except ValueError:
        print("Введенные данные не преобразуются в строку")
    try:
        if len(inp_str1) != len(inp_str2) or inp_str.shape != tuple([2]):
            raise Exception()
        for i in range(len(inp_str1)):
            if DNA_compare.get(inp_str1[i]) != inp_str2[i] or DNA_compare.get(inp_str1[i]) is None:
                raise Exception()
    except ValueError:
        print("Введенные данные не являются ДНК")
    return np.array([inp_str1, inp_str2])


def is_rna(inp_str):
    try:
        inp_str = str(inp_str)
    except ValueError:
        print("Введенные данные не преобразуются в строку")
    try:
        for i in range(len(inp_str)):
            if RNA_compare.get(inp_str[i]) is None:
                raise Exception()
    except ValueError:
        print("Введенные данные не являются РНК")
    return inp_str


def rna_transcription(inp_rna):
    res = ""
    for i in range(len(inp_rna)):
        res += RNA_compare[inp_rna[i]]
    return res


def complimentary_dna(inp_dna):
    res = ""
    for i in range(len(inp_dna)):
        res += DNA_compare[inp_dna[i]]
    return res


def mult(f_inp, s_inp):
    rnd_numbers = np.random.randint(2, size=min(len(f_inp), len(s_inp)))
    tmp_res = ""
    for i in range(len(rnd_numbers)):
        if rnd_numbers[i] == 0:
            tmp_res += f_inp[i]
        else:
            tmp_res += s_inp[i]
    if len(f_inp) > len(s_inp):
        tmp_res += f_inp[len(s_inp):]
    else:
        tmp_res += s_inp[len(f_inp):]
    return tmp_res


class sequence(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __getitem__(self, item):
        pass

    @abstractmethod
    def __setitem__(self, key, value):
        pass

    @abstractmethod
    def __add__(self, other):
        pass

    @abstractmethod
    def __mul__(self, other):
        pass

    @abstractmethod
    def __eq__(self, other):
        pass

    @abstractmethod
    def __str__(self):
        pass


class DNA(sequence, ABC):
    def __init__(self, inp_str):
        self.array = is_dna(inp_str)

    def __getitem__(self, item):
        try:
            return np.array([self.array[0][item], self.array[1][item]])
        except ValueError:
            print("Длина ДНК меньше введенного числа")

    def get(self):
        return np.array([self.array[0], self.array[1]])

    def __setitem__(self, key, value):
        tmp = is_dna(value)
        try:
            if len(tmp[0]) != 1:
                raise Exception()
        except ValueError:
            print("Введенные данные не являются одним элементом ДНК")
        try:
            self.array[0][key] = is_dna(value)[0][0]
            self.array[1][key] = is_dna(value)[1][0]
        except ValueError:
            print("Длина ДНК меньше введенного числа")

    def __add__(self, other):
        try:
            if isinstance(other, DNA):
                tmp = other.get()
                return DNA(np.array([self.array[0] + tmp[0], self.array[1] + tmp[1]]))
            else:
                raise Exception()
        except ValueError:
            print("Объекты разных типов данных")

    def __mul__(self, other):
        try:
            if isinstance(other, DNA):
                tmp = other.get()
                tmp_m = mult(self.array[0], tmp[0])
                return DNA(np.array([tmp_m, complimentary_dna(tmp_m)]))
            else:
                raise Exception()
        except ValueError:
            print("Объекты разных типов данных")

        other = is_dna(other)
        first_sequence = mult(self.array[0], other)
        self.array[0] = first_sequence
        self.array[1] = complimentary_dna(first_sequence)

    def __eq__(self, other):
        if isinstance(other, DNA):
            return self.array == other.get()
        return False

    def __str__(self):
        return np.array2string(self.array, separator=' ')


class RNA(sequence, ABC):
    def __init__(self, inp_str):
        self.array = is_rna(inp_str)

    def __getitem__(self, item):
        try:
            return self.array[item]
        except ValueError:
            print("Длина РНК меньше введенного числа")

    def __setitem__(self, key, value):
        tmp = is_rna(value)
        try:
            if len(tmp) != 1:
                raise Exception()
            self.array[key] = is_rna(value)
        except ValueError:
            print("Введенные данные не являются одним элементом РНК")

    def __add__(self, other):
        try:
            if isinstance(other, RNA):
                return RNA(self.array + str(other))
            else:
                raise Exception()
        except ValueError:
            print("Объекты разных типов данных")

    def __mul__(self, other):
        try:
            if isinstance(other, RNA):
                return RNA(mult(self.array, str(other)))
            else:
                raise Exception()
        except ValueError:
            print("Объекты разных типов данных")

    def __eq__(self, other):
        if isinstance(other, RNA):
            return self.array == str(other)
        return False

    def __str__(self):
        return self.array

    def build_dna(self):
        tmp = rna_transcription(self.array)
        return DNA(np.array([tmp, complimentary_dna(tmp)]))


if __name__ == '__main__':
    a = RNA("AUUGAACUA")
    b = RNA("CGGAAA")
    print(a, b)
    print(a[2], b[4])
    print(a + b)
    print(a * b)
    A = a.build_dna()
    B = b.build_dna()
    print(A, B)
    print(A[4], B[0])
    print(A + B)
    print(A * B)

a = RNA("AUUGAACUA")
b = RNA("CGGAAA")
print(a + b)