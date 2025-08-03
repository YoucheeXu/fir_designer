#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import sys
import re
import inspect
import time
from enum import IntEnum
from functools import wraps, partial
from typing import override

import numpy as np


def pv3(var, var_name: str):
    """ print 变量的名字及type、及其值 """
    print(f"{var_name} : {type(var).__name__}; value: ")
    print(var)


# namestr(a, globals())
def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]


def ps(var, var_name: str):
    """ print 变量的名字及type、shape """
    print(f"{var_name} : {var.__class__.__name__}; shape: ")
    print(np.shape(var))


# 这是从堆栈内存的第3层开始查找返回变量名称
def retrieve_name(var):
    # print(f"globals = <{globals()}>")
    # print(f"locals = <{locals()}>")
    print(f'var = {var}')
    for fi in inspect.stack()[2:]:
        for item in fi.frame.f_locals.items():
            if var is item[1]:
                return item[0]
    return ""


def get_variable_name(variable):
    loc = locals()
    for key in loc:
        if loc[key] == variable:
            return key


#### untest above ####


def fcn_name(fcn):
    str_list = str(fcn).split()
    # print(str_list)
    return str_list[1]


def po(*values: object, endstr: str = "\n"):
    f_back = sys._getframe().f_back
    filename = f_back.f_code.co_filename
    linenum = f_back.f_lineno

    print(f"{linenum}@{filename} ", end="")
    for val in values:
        print(val, end=", ")
    print(end=endstr)

# 1 dimision array: ary[i]                  OK
# 2 dimision array: ary[i, j] ary[i][j]     OK
# slice
def pv(p: object, endstr: str = "\n"):
    var_name = ""
    for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
        if (m := re.search(r'\bpv\s*\((.+)\)', line)):
            var_name = m.group(1).split(", end")[0]
            break
    i_reg = r"([^,:]+)"
    interval = ""
    if(m := re.search(rf"(.+\[){i_reg}\]\[{i_reg}\]", var_name)):
        interval = "]["
    elif(m := re.search(rf"(.+\[){i_reg},{i_reg}\]", var_name)):
        interval = ", "
    elif(m := re.search(rf"(.+\[){i_reg}\]", var_name)):
        interval = ""
    else:
        pass

    if m:
        try:
            item_class = inspect.currentframe().f_back.f_locals.items()
            # print(itemClass)
            mgroups = m.groups()
            # print(mGroups)
            # varName = m.group(1)
            var_name = mgroups[0]
            len_mgroups = len(mgroups)
            idx_lst = [mgroups[i].strip() for i in range(1, len_mgroups)]
            # print(idx_lst)
            idx_dict = {var_name: str(var_val) for var_name, var_val in item_class if var_name in idx_lst}
            # print(idx_dict)
            val_lst = [idx_dict.get(idx_lst[i - 1], idx_lst[i - 1]) for i in range(1, len_mgroups)]
            var_name += interval.join(val_lst) + "]"
        except Exception as r:
            print("error =", r)

    f_back = sys._getframe().f_back
    filename = f_back.f_code.co_filename
    linenum = f_back.f_lineno

    print(f'{linenum}@{filename} {var_name} = {p}', end=endstr)


def pe(exp, end: str = "\n"):
    # exp_name = str(exp)
    # i = 0
    for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
        # print(i, line)
        if m := re.search(r"\bpe\s*\((.+)\)", line):
            # print(m.groups())
            exp_name = m.group(1)
            # var = m.group(2)
            # print(f'{exp} = {var}', end=end)
            break
        # i += 1
    if exp:
        pass
    print(f'{exp_name} = {exp}', end=end)


# pe = partial(pe2, var="")

# 定义装饰器
def time_calc(func):
    @wraps(func)
    def wrapper(*args, **kargs):
        start_time = time.time()
        f = func(*args, **kargs)
        exec_time = time.time() - start_time
        print(f"{func.__name__}花费的时间是：{exec_time}")
        return f
    return wrapper


class LogLevel(IntEnum):
    INFO = 1
    WARN = 2
    ERROR = 3


class Logit():
    def __init__(self, level: LogLevel = LogLevel.INFO, logfile: str = ""):
        self.__level: LogLevel = level
        self.__logfile: str = logfile

    def __call__(self, func): # 接受函数
        @wraps(func)
        def wrapper(*args, **kwargs):
            # print(__name__)
            # fun_name = sys._getframe().f_code.co_name
            fun_name = func.__name__
            self._log(self.__level.name, f"{fun_name}() was called")
            return func(*args, **kwargs)
        return wrapper  # 返回函数

    def _notify(self):
        pass

    def _log(self, level: str, msg: str):
        # print(__file__)
        # print(sys._getframe().f_lineno)
        timestr = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
        filename = sys._getframe().f_code.co_filename
        linenum = sys._getframe().f_back.f_back.f_lineno
        logstring = f"{timestr} {filename}@{linenum:03d} [{level}]: {msg}"
        print(logstring)
        if self.__logfile:
            with open(self.__logfile, 'a', encoding='utf8') as opened_file:
                _ = opened_file.write(logstring + '\n')
        self._notify()

    def log_info(self, msg: str):
        self._log("INFO", msg)

    def log_warn(self, msg: str):
        self._log("Warning", msg)

    def log_err(self, msg: str):
        self._log("Error", msg)


class EmailLogit(Logit):
    def __init__(self, email: str, *args, **kwargs):
        self.__email = email
        super().__init__(*args, **kwargs)

    @override
    def _notify(self):
        # send a email to self.__email
        pass


if __name__ == '__main__':

    # import os
    # os.system('pause')

    # @time_calc
    # def want_sleep(sleep_time: float):
        # time.sleep(sleep_time)

    # @Logit()
    # def say(something: str):
        # print(f"say {something}!")

    # want_sleep(1)
    # say('hello')
    # logit = Logit()
    # logit.log_info("haha")
    # logit.log_warn("hehe")

    '''
    # a = 1

    # val = namestr(a, globals())
    # print(val)
    # print(retrieve_name(a))

    # pv3(a + 1, a)

    # pv(a)
    # pv(a + 2, end="\t")
    # pv(a + 3)

    # import numpy as np
    # ary1 = np.linspace(0, 39, 40)

    # ary2 = [i for i in range(40)]

    # for i in range(40):
        # pv(ary1[i])
        # pv(ary2[i])

    # i = 1

    # pv(ary1[i])
    # pv(ary2[i])

    # i = 1
    # j = 2
    # matij = np.random.randn(3, 3)
    # pv(matij[i, j])
    # pv(matij[i][j])


    # str1 = "Hello World"
    # m = re.findall('(Hello)|(World)', str1)
    # print(m)
    # print(m.groups())
    '''
    pe(2 / 5)
