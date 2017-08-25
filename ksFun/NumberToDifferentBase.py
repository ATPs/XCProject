# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 13:30:12 2015

@author: k
given a number and a base, the default is 36, return a string based on the base.
for example, 10: A, 11:B, 12:C
"""

def baseEncode(number, base=36):
    """
    given a int number, return a string based on base
    10:A, 36:Z, 37:Z1, under default setting
    """
    if not isinstance(number, int):
        raise TypeError('number must be an integer')
    alphabet='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    if base > 62 or base <=1:
        print("base should be between 2 and 62")
        return None
    sign = ""
    if number < 0:
        sign = "-"
        number = -number
    alphabet = alphabet[:base+1]
    if 0 <= number and number <base:
        return sign+alphabet[number]
    numberbase=""
    while number != 0:
        number, i = divmod(number, base)
        numberbase = alphabet[i] + numberbase
    return sign+numberbase

def reverseBaseEncode(numberstr, base=36):
    """
    given a string of baseEncode, return a int of it represent
    """
    if not isinstance(numberstr, str):
        raise TypeError('number must be an str')
    alphabet='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    alphabetDic = {}
    for i in range(base):
        alphabetDic[alphabet[i]] = i
        
    if base > 62 or base <=1:
        print("base should be between 2 and 62")
        return None
    sign = "+"
    if numberstr[0] =="-" or numberstr[0] == "+":
        sign = numberstr[0]
        numberstr = numberstr[1:]
    number = alphabetDic[numberstr[0]]
    for ele in numberstr[1:]:
        number = number*base + alphabetDic[ele]
    if sign == "-":
        number = -number
    return number
    