#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import sys
import os
import math
from enum import IntEnum, unique, Enum

import tkinter as tk
import tkinter.messagebox as tkMessageBox
import tkinter.filedialog as tkFileDialog
import tkinter.ttk as ttk

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk #NavigationToolbar2TkAgg
from matplotlib.pylab import mpl

from logit import pv


np.seterr(divide='ignore', invalid='ignore')


@unique
class WinType(IntEnum):
	Rectangular = 0
	Triangular = 1
	Welch = 2
	Sine = 3
	Hann = 4
	Hamming = 5
	Blackman = 6
	Nuttall = 7
	BlackmanNuttall = 8
	BlackmanHarris = 9
	FlatTop = 10

@unique
class FilrType(Enum):
	LowPass = 0
	HighPass = 1
	BandPass = 2
	BandStop = 3

@unique
class RespPltType(Enum):
	Impulse = 0
	Step = 1


class MatPlotData:
    def __init__(self):
        self.canvas = None
        self.ax = None
        self.lines = None

class App(tk.Frame):

	def __init__(self, master):
		super().__init__(master)

		# Initial settings

		# Parameters
		self.__tiSample = 0.01
		self.__tiSampleVar = tk.DoubleVar()
		self.__tiSampleVar.set(self.__tiSample)

		self.__frqCutOffLo = 10.0
		self.__frqCutOffLoVar = tk.DoubleVar()
		self.__frqCutOffLoVar.set(self.__frqCutOffLo)

		self.__frqCutOffHi = 20.0
		self.__frqCutOffHiVar = tk.DoubleVar()
		self.__frqCutOffHiVar.set(self.__frqCutOffHi)

		self.__numTotalSample = 64
		self.__numTotalSampleVar = tk.IntVar()
		self.__numTotalSampleVar.set(self.__numTotalSample)

		self.__numShftSample = 32
		self.__numShftSampleVar = tk.IntVar()
		self.__numShftSampleVar.set(self.__numShftSample)

		self.__stWinTyp = "Hamming"
		self.__stWinTypVar = tk.StringVar()
		self.__stWinTypVar.set(self.__stWinTyp)

		self.__stFiltTyp = FilrType.LowPass
		self.__stFiltTypVar = tk.IntVar()
		self.__stFiltTypVar.set(0)	# FilrType.LowPass

		self.__stRespPltTyp = RespPltType.Impulse
		self.__stRespPltTypVar = tk.IntVar()
		self.__stRespPltTypVar.set(0)	# Impulse

		self.__stDisplay = tk.StringVar()
		self.__stDisplayVar = tk.StringVar()
		self.__stDisplayVar.set("Both")

		self.__numFreqSample = 256

		# Time Domain
		numTotalSample = self.__numTotalSample
		self.__timeVec = [0.0] * self.__numTotalSample
		self.__impulseResp = [0.0] * self.__numTotalSample
		self.__stepResp = [0.0] * self.__numTotalSample
		self.__window = [0.0] * self.__numTotalSample
		self.__windowedImpResp = [0.0] * self.__numTotalSample
		self.__windowedStepResp = [0.0] * self.__numTotalSample

		# Frequency Domain
		numFreqSample = self.__numFreqSample
		self.__freqVec = [0.0] * numFreqSample
		self.__impRespMag = [0.0] * numFreqSample
		self.__winRespMag = [0.0] * numFreqSample
		self.__winMag = [0.0] * numFreqSample

		self.__CreateWindow(master)

		# cmbWin.SelectedIndex = 5
		# radViewImpulse.Checked = True
		# cmbDisplay.SelectedIndex = 0

		self.__FilrTypeChange()

		self.__DesignFilter()

	def __ShowErr(self, err):
		tkMessageBox.showerror("Design FIR", err)

	def __ShowInfo(self, title, info):
		tkMessageBox.showinfo(title, info)

	'''
		设置窗口居中和宽高
		:param window: 主窗体
		:param Width: 窗口宽度
		:param Hight: 窗口高度
		:return: 无
	'''
	def __CenterWindow(self, window, width, hight):

		# 获取屏幕宽度和高度
		sw = window.winfo_screenwidth()
		sh = window.winfo_screenheight()

		# 计算中心坐标
		cen_x = (sw - width) / 2
		cen_y = (sh - hight) / 2 * 0.9

		# 设置窗口大小并居中
		window.geometry('%dx%d+%d+%d' %(width, hight, cen_x, cen_y))

	def __CreateWindow(self, root):

		root.title("Design FIR")  # 设置窗口标题
		self.__CenterWindow(root, 1200, 500)

		frmTiDomain = tk.Frame(root)

		# Filter Response in Time Domain
		self.__dataFilrTiDomain = MatPlotData()
		figFilrTiDomain, self.__dataFilrTiDomain.ax, self.__dataFilrTiDomain.lines = self.__CreateMatplot(r"Filter Response in Time Domain", r"Time (s)", r"", 4)
		self.__dataFilrTiDomain.canvas = FigureCanvasTkAgg(figFilrTiDomain, frmTiDomain)
		self.__dataFilrTiDomain.canvas.draw()
		self.__dataFilrTiDomain.canvas.get_tk_widget().pack(side=tk.LEFT, fill='both', expand=True, padx=(10,10), pady=5)
		# canvasFilrTiDomain.get_tk_widget().grid(row=0,column=0, columnspan = 4)

		# Window Function in Time Domain
		self.__dataWinTiDomain = MatPlotData()
		figWindowTiDomain, self.__dataWinTiDomain.ax, self.__dataWinTiDomain.lines = self.__CreateMatplot(r"Window Function in Time Domain", r"Time (s)", r"")
		self.__dataWinTiDomain.canvas = FigureCanvasTkAgg(figWindowTiDomain, frmTiDomain)
		self.__dataWinTiDomain.canvas.draw()
		self.__dataWinTiDomain.canvas.get_tk_widget().pack(side=tk.LEFT, fill='both', expand=True, padx=(10,10), pady=5)
		# canvasWindowTiDomain.get_tk_widget().grid(row=0, column=4, columnspan = 6)

		frmTiDomain.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=(10,5))

		frmFreqDomain = tk.Frame(root)

		# Filter Frequency Response
		self.__dataFilrFreqDomain = MatPlotData()
		figFilrFreqDomain, self.__dataFilrFreqDomain.ax, self.__dataFilrFreqDomain.lines = self.__CreateMatplot(r"Filter Frequency Response", r"Frequency (Hz)", r"", 2)
		self.__dataFilrFreqDomain.canvas = FigureCanvasTkAgg(figFilrFreqDomain, frmFreqDomain)
		self.__dataFilrFreqDomain.canvas.draw()
		self.__dataFilrFreqDomain.canvas.get_tk_widget().pack(side=tk.LEFT, fill='both', expand=True, padx=(10,10), pady=5)
		# canvasFilrFrqDomain.get_tk_widget().grid(row=1,column=0, columnspan = 4)

		# Window Function Frequency Response
		self.__dataWinFreqDomain = MatPlotData()
		figWindowFreqDomain, self.__dataWinFreqDomain.ax, self.__dataWinFreqDomain.lines = self.__CreateMatplot(r"Window Function Frequency Response", r"Frequency (Hz)", r"")
		self.__dataWinFreqDomain.canvas = FigureCanvasTkAgg(figWindowFreqDomain, frmFreqDomain)
		self.__dataWinFreqDomain.canvas.draw()
		self.__dataWinFreqDomain.canvas.get_tk_widget().pack(side=tk.LEFT, fill='both', expand=True, padx=(10,10), pady=5)
		# canvasWindowFrqDomain.get_tk_widget().grid(row=1, column=4, columnspan = 6)

		frmFreqDomain.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=(5,5))

		frmParas = tk.LabelFrame(root, text="Parameters")

		lblSampTi = tk.Label(frmParas, text = "Sample Time(s):")
		lblSampTi.grid(row=0,column=0)
		enySampTi = ttk.Entry(frmParas, textvariable = self.__tiSampleVar, width = 15)
		enySampTi.grid(row=0,column=1)

		lblFilrLength = tk.Label(frmParas, text = "Filter Length:")
		lblFilrLength.grid(row=1,column=0)
		enyFilrLength = ttk.Entry(frmParas, textvariable = self.__numTotalSampleVar, width = 15)
		enyFilrLength.grid(row=1,column=1)

		lblShiftSamples = tk.Label(frmParas, text = "Shift Samples:")
		lblShiftSamples.grid(row=2,column=0)
		enyShiftSamples = ttk.Entry(frmParas, textvariable = self.__numShftSampleVar, width = 15)
		enyShiftSamples.grid(row=2,column=1)

		lblCutOffFrqLo = tk.Label(frmParas, text = "Low Frequency(Hz):")
		lblCutOffFrqLo.grid(row=0,column=2)
		self.__txtCutOffFrqLo = ttk.Entry(frmParas, textvariable = self.__frqCutOffLoVar, width = 15)
		self.__txtCutOffFrqLo.grid(row=0,column=3)

		lblCutOffFrqHi = tk.Label(frmParas, text = "High Frequency(Hz):")
		lblCutOffFrqHi.grid(row=1,column=2)
		self.__txtCutOffFrqHi = ttk.Entry(frmParas, textvariable = self.__frqCutOffHiVar, width = 15)
		self.__txtCutOffFrqHi.grid(row=1,column=3)

		lblcmbWinTyp = tk.Label(frmParas, text = "Window Function:")
		lblcmbWinTyp.grid(row=2,column=2)
		options = ["Rectangular", "Triangular", "Welch", "Sine", "Hann", "Hamming", "Blackman", "Nuttall", "BlackmanNuttall", "BlackmanHarris", "FlatTop"]
		cmbWinTyp = ttk.Combobox(frmParas, textvariable = self.__stWinTypVar, values = options, state='readonly', width = 15)
		cmbWinTyp.current(5)
		cmbWinTyp.bind('<<ComboboxSelected>>', self.__cmbWin_SelectedIndexChanged)
		cmbWinTyp.grid(row=2,column=3)

		# radBS radBP radHP radLP
		lblFiltTyp = tk.Label(frmParas, text = "Filter Type")
		lblFiltTyp.grid(row=0,column=4, columnspan = 2)
		radio = ttk.Radiobutton(frmParas, variable = self.__stFiltTypVar, value = 0, text ='Low Pass', command = self.__radPS_CheckedChanged)
		radio.grid(row = 1,column = 4, sticky = tk.W)
		radio = ttk.Radiobutton(frmParas, variable = self.__stFiltTypVar, value = 1, text ='High Pass', command = self.__radPS_CheckedChanged)
		radio.grid(row=1,column=5, sticky = tk.W)
		radio = ttk.Radiobutton(frmParas, variable = self.__stFiltTypVar, value = 2, text ='Band Pass', command = self.__radPS_CheckedChanged)
		radio.grid(row=2,column=4, sticky = tk.W)
		radio = ttk.Radiobutton(frmParas, variable = self.__stFiltTypVar, value = 3, text ='Band Stop', command = self.__radPS_CheckedChanged)
		radio.grid(row=2,column=5, sticky = tk.W)

		# radViewStep radViewImpulse
		lblRespPltTyp = tk.Label(frmParas, text = "Response Plot:")
		lblRespPltTyp.grid(row=0,column=6)
		radio = ttk.Radiobutton(frmParas, variable = self.__stRespPltTypVar, value = 0, text ='Impulse', command = self.__radViewImpulse_CheckedChanged)
		radio.grid(row=0,column=7, sticky = tk.W)
		radio = ttk.Radiobutton(frmParas, variable = self.__stRespPltTypVar, value = 1, text ='Step', command = self.__radViewImpulse_CheckedChanged)
		radio.grid(row=1,column=7, sticky = tk.W)

		# cmbDisplay
		lblcmbDisplay = tk.Label(frmParas, text = "Display:")
		lblcmbDisplay.grid(row=2,column=6)
		options = ["Both", "Raw Only", "Win Only"]
		cmbDisplay = ttk.Combobox(frmParas, textvariable = self.__stDisplayVar, values = options, state='readonly', width = 15)
		cmbDisplay.current(0)
		cmbDisplay.bind('<<ComboboxSelected>>', self.__cmbDisplay_SelectedIndexChanged)
		cmbDisplay.grid(row=2,column=7)

		frmParas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(10,5), pady=(5,10))

		btnDesignFilr = ttk.Button(master = root, text = "Design Filter", command = self.__btnDesignFilr_Click)
		# btnDesignFilr.grid(row = 2,column = 8, rowspan = 3, columnspan = 2, sticky = tk.N + tk.S + tk.W + tk.E)
		btnDesignFilr.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(5,10), pady=10)


		# 创建一个顶级菜单
		menubar = tk.Menu(root)

		# 在顶层窗口menubar中加入下拉窗口filemenu
		filemenu = tk.Menu(menubar, tearoff=False)

		filemenu.add_command(label="Export Coefficients...", command=self.__exportCoefficientsToolStripMenuItem_Click)
		filemenu.add_command(label="Export Time Domain Data...", command=self.__exportTiDomainDataToolStripMenuItem_Click)
		filemenu.add_command(label="Export Frequency Domain Data...", command=self.__exportFreqDomainDataToolStripMenuItem_Click)
		filemenu.add_separator()

		filemenu.add_command(label="Exit", command=root.destroy)
		menubar.add_cascade(label="File", menu=filemenu)	# 为filemenu命名为‘文件’

		# infoToolStripMenuItem
		menubar.add_command(label = "Info", command=self.__infoToolStripMenuItem_Click)

		# 显示菜单
		root.config(menu=menubar)

	def __CreateMatplot(self, title, xLabel, yLabel, numLine = 1):

		fig = plt.Figure(figsize=(5.8, 1.8), dpi=100)
		ax = fig.add_subplot()
		fig.subplots_adjust(bottom=0.25)

		t = np.arange(0, 3, 0.01)
		y = np.sin(2 * np.pi * t)

		line = []
		for n in range(numLine):
			l, = ax.plot(t, y)
			line.append(l)

		ax.set_title(title)
		ax.set_xlabel(xLabel)
		ax.set_ylabel(yLabel)
		ax.grid()
		ax.autoscale()

		return fig, ax, line

	def __DesignFilter(self):
		self.__ComputeTimeVec()
		self.__ComputeWindow()
		self.__ComputeResps()
		self.__ComputeWindowedResps()

		self.__ComputeFreqVec()
		self.__ComputeRespBode()
		self.__ComputeWindowDFT()

		self.__UpdateCharts()
		self.__UpdatePlotSettings()

	def __UpdateCharts(self):

		'''
		chrtFilrTiDomain.lines[0].Points.DataBindXY(self.__timeVec, self.__impulseResp)
		chrtFilrTiDomain.lines[1].Points.DataBindXY(self.__timeVec, self.__windowedImpResp)
		chrtFilrTiDomain.lines[2].Points.DataBindXY(self.__timeVec, self.__stepResp)
		chrtFilrTiDomain.lines[3].Points.DataBindXY(self.__timeVec, self.__windowedStepResp)
		chrtFilrTiDomain.lines[2].set_visible(False)
		chrtFilrTiDomain.lines[3].set_visible(False)

		chrtFilrTiDomain.ChartAreas[0].RecalculateAxesScale()
		chrtFilrTiDomain.Update()
		'''
		self.__dataFilrTiDomain.lines[0].set_data(self.__timeVec, self.__impulseResp)
		self.__dataFilrTiDomain.lines[1].set_data(self.__timeVec, self.__windowedImpResp)
		self.__dataFilrTiDomain.lines[2].set_data(self.__timeVec, self.__stepResp)
		self.__dataFilrTiDomain.lines[2].set_visible(False)
		self.__dataFilrTiDomain.lines[3].set_data(self.__timeVec, self.__windowedStepResp)
		self.__dataFilrTiDomain.lines[3].set_visible(False)
		self.__dataFilrTiDomain.ax.set_xlim(self.__timeVec[0], self.__timeVec[-1])
		yMin = min(min(self.__impulseResp), min(self.__windowedImpResp))
		yMax = max(max(self.__impulseResp), max(self.__windowedImpResp))
		self.__dataFilrTiDomain.ax.set_ylim(yMin, yMax)
		self.__dataFilrTiDomain.canvas.draw()

		'''
		chrtWinTiDomain.lines[0].Points.DataBindXY(self.__timeVec, self.__window)
		chrtWinTiDomain.ChartAreas[0].RecalculateAxesScale()
		'''
		self.__dataWinTiDomain.lines[0].set_data(self.__timeVec, self.__window)
		self.__dataWinTiDomain.ax.set_xlim(self.__timeVec[0], self.__timeVec[-1])
		self.__dataWinTiDomain.ax.set_ylim(min(self.__window), max(self.__window))
		self.__dataWinTiDomain.canvas.draw()

		'''
		chrtFilrFreqDomain.lines[0].Points.DataBindXY(self.__freqVec, self.__impRespMag)
		chrtFilrFreqDomain.lines[1].Points.DataBindXY(self.__freqVec, self.__winRespMag)
		chrtFilrFreqDomain.ChartAreas[0].AxisX.Interval = 10.0 * (0.5 / self.__tiSample) / (self.__numFreqSample - 1.0)
		chrtFilrFreqDomain.ChartAreas[0].RecalculateAxesScale()
		'''

		self.__dataFilrFreqDomain.lines[0].set_data(self.__freqVec, self.__impRespMag)
		self.__dataFilrFreqDomain.lines[1].set_data(self.__freqVec, self.__winRespMag)
		self.__dataFilrFreqDomain.ax.set_xlim(self.__freqVec[0], self.__freqVec[-1])
		yMin = min(min(self.__impRespMag), min(self.__winRespMag))
		yMax = max(max(self.__impRespMag), max(self.__winRespMag))
		self.__dataFilrFreqDomain.ax.set_ylim(yMin, yMax)
		self.__dataFilrFreqDomain.canvas.draw()

		'''
		chrtWinFreqDomain.lines[0].Points.DataBindXY(self.__freqVec, self.__winMag)
		chrtWinFreqDomain.ChartAreas[0].AxisX.Interval = 10.0 * (0.5 / self.__tiSample) / (self.__numFreqSample - 1.0)
		chrtWinFreqDomain.ChartAreas[0].RecalculateAxesScale()
		'''
		self.__dataWinFreqDomain.lines[0].set_data(self.__freqVec, self.__winMag)
		self.__dataWinFreqDomain.ax.set_xlim(self.__freqVec[0], self.__freqVec[-1])
		self.__dataWinFreqDomain.ax.set_ylim(min(self.__winMag), max(self.__winMag))
		self.__dataWinFreqDomain.canvas.draw()


	# Time Domain Functions
	def __ComputeTimeVec(self):
		numTotalSample = self.__numTotalSample
		# self.__timeVec = new float[self.__numTotalSample]
		self.__timeVec = [0.0] * numTotalSample

		tiSample = self.__tiSample
		# for (int n = 0 n < self.__numTotalSample n++):
		for n in range(numTotalSample):
			self.__timeVec[n] = n * tiSample

		# pv(self.__timeVec)

	def __ComputeResps(self):
		numTotalSample = self.__numTotalSample
		# self.__impulseResp = new float[self.__numTotalSample]
		self.__impulseResp = [0.0] * self.__numTotalSample
		# self.__stepResp = new float[self.__numTotalSample]
		self.__stepResp = [0.0] * self.__numTotalSample

		numShftSample = self.__numShftSample
		stFiltTyp = self.__stFiltTyp
		tiSample = self.__tiSample
		frqCutOffLo = self.__frqCutOffLo
		frqCutOffHi = self.__frqCutOffHi
		# for (int n = 0 n < self.__numTotalSample n++)
		for n in range(numTotalSample):
			if (n != numShftSample):
				if stFiltTyp == FilrType.LowPass:
					self.__impulseResp[n] = math.sin(2.0 * math.pi * frqCutOffLo * tiSample * (n - numShftSample)) / (math.pi * tiSample * (n - numShftSample))
				elif stFiltTyp == FilrType.HighPass:
					self.__impulseResp[n] = (math.sin(math.pi * (n - numShftSample)) - math.sin(2.0 * math.pi * frqCutOffHi * tiSample * (n - numShftSample))) / (math.pi * tiSample * (n - numShftSample))
				elif stFiltTyp == FilrType.BandPass:
					self.__impulseResp[n] = (math.sin(2.0 * math.pi * frqCutOffHi * tiSample * (n - numShftSample)) - math.sin(2.0 * math.pi * frqCutOffLo * tiSample * (n - numShftSample))) / (math.pi * tiSample * (n - numShftSample))
				elif stFiltTyp == FilrType.BandStop:
					self.__impulseResp[n] = (math.sin(2.0 * math.pi * frqCutOffLo * tiSample * (n - numShftSample)) - math.sin(2.0 * math.pi * frqCutOffHi * tiSample * (n - numShftSample)) + math.sin(math.pi * (n - numShftSample))) / (math.pi * tiSample * (n - numShftSample))

			else: # Avoid divide-by-zero, limit is 2*fc
				if stFiltTyp == FilrType.LowPass:
					self.__impulseResp[n] = 2.0 * frqCutOffLo
				elif stFiltTyp == FilrType.HighPass:
					self.__impulseResp[n] = 1.0 / tiSample - 2.0 * frqCutOffHi
				elif stFiltTyp == FilrType.BandPass:
					self.__impulseResp[n] = 2.0 * frqCutOffHi - 2.0 * frqCutOffLo
				elif stFiltTyp == FilrType.BandStop:
					self.__impulseResp[n] = 2.0 * frqCutOffLo - 2.0 * frqCutOffHi + 1.0 / tiSample

		# Normalise by DC gain to achieve 0dB gain at DC and then compute step response
		# for (int n = 0 n < self.__numTotalSample n++)
		for n in range(numTotalSample):
			self.__impulseResp[n] *= tiSample

			if (n == 0):
				self.__stepResp[n] = 0.5 * self.__impulseResp[n]
			else:
				self.__stepResp[n] = self.__stepResp[n - 1] + 0.5 * (self.__impulseResp[n] + self.__impulseResp[n - 1])

		# pv(self.__impulseResp)
		# pv(self.__stepResp)

	def __ComputeWindow(self):
		numTotalSample = self.__numTotalSample
		# self.__window = new float[self.__numTotalSample]
		self.__window = [0.0] * self.__numTotalSample

		# for (int n = 0 n < self.__numTotalSample n++)
		for n in range(numTotalSample):
			if self.__stWinTyp == "Rectangular":
				self.__window[n] = 1.0
			elif self.__stWinTyp == "Triangular":
				self.__window[n] = 1.0 - abs((n - 0.5 * self.__numTotalSample) / (0.5 * self.__numTotalSample))
			elif self.__stWinTyp == "Welch":
				self.__window[n] = 1.0 - math.pow((n - 0.5 * self.__numTotalSample) / (0.5 * self.__numTotalSample), 2.0)
			elif self.__stWinTyp == "Sine":
				self.__window[n] = math.sin(math.pi * n / (numTotalSample))
			elif self.__stWinTyp == "Hann":
				self.__window[n] = 0.5 * (1 - math.cos(2.0 * math.pi * n / (numTotalSample)))
			elif self.__stWinTyp == "Hamming":
				self.__window[n] = (25.0 / 46.0) - (21.0 / 46.0) * math.cos(2.0 * math.pi * n / (numTotalSample))
			elif self.__stWinTyp == "Blackman":
				self.__window[n] = 0.42 - 0.5 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.08 * math.cos(4.0 * math.pi * n / (numTotalSample))
			elif self.__stWinTyp == "Nuttall":
				self.__window[n] = 0.355768 - 0.487396 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.144232 * math.cos(4.0 * math.pi * n / (numTotalSample)) - 0.012604 * math.cos(6.0 * math.pi * n / (numTotalSample))
			elif self.__stWinTyp == "BlackmanNuttall":
				self.__window[n] = 0.3635819 - 0.4891775 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.1365995 * math.cos(4.0 * math.pi * n / (numTotalSample)) - 0.0106411 * math.cos(6.0 * math.pi * n / (numTotalSample))
			elif self.__stWinTyp == "BlackmanHarris":
				self.__window[n] = 0.35875 - 0.48829 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.14128 * math.cos(4.0 * math.pi * n / (numTotalSample)) - 0.01168 * math.cos(6.0 * math.pi * n / (numTotalSample))
			elif self.__stWinTyp == "FlatTop":
				self.__window[n] = 0.21557895 - 0.41663158 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.277263158 * math.cos(4.0 * math.pi * n / (numTotalSample)) - 0.083578947 * math.cos(6.0 * math.pi * n / (numTotalSample)) + 0.006947368 * math.cos(8.0 * math.pi * n / (numTotalSample))
			else:
				self.__window[n] = 1.0

		# pv(self.__window)

	def __ComputeWindowedResps(self):
		numTotalSample = self.__numTotalSample
		# self.__windowedImpResp = new float[self.__numTotalSample]
		self.__windowedImpResp = [0.0] * self.__numTotalSample
		# self.__windowedStepResp = new float[self.__numTotalSample]
		self.__windowedStepResp = [0.0] * self.__numTotalSample

		# for (int n = 0; n < self.__numTotalSample; n++):
		for n in range(numTotalSample):
			self.__windowedImpResp[n] = self.__impulseResp[n] * self.__window[n]

			if (n == 0):
				self.__windowedStepResp[n] = 0.5 * self.__windowedStepResp[n]
			else:
				self.__windowedStepResp[n] = self.__windowedStepResp[n - 1] + 0.5 * (self.__windowedImpResp[n] + self.__windowedImpResp[n - 1])

		# pv(self.__windowedImpResp)
		# pv(self.__windowedStepResp)

	# Frequency Domain Functions
	def __ComputeFreqVec(self):
		numFreqSample = self.__numFreqSample
		# freqVecHz = new float[self.__numFreqSample]
		self.__freqVec = [0.0] * numFreqSample

		tiSample = self.__tiSample
		df = (0.5 / tiSample) / (numFreqSample - 1.0)
		pv(df)

		# for (int n = 0; n < self.__numFreqSample; n++):
		for n in range(numFreqSample):
			self.__freqVec[n] = n * df

	def __ComputeRespBode(self):
		numFreqSample = self.__numFreqSample
		# impRespMag = new float[self.__numFreqSample]
		self.__impRespMag = [0.0] * numFreqSample
		# winRespMag = new float[self.__numFreqSample]
		self.__winRespMag = [0.0] * numFreqSample

		numTotalSample = self.__numTotalSample
		tiSample = self.__tiSample
		# for (int fIndex = 0; fIndex < NUM_FREQ_SAMPLES; fIndex++):
		for fIndex in range(numFreqSample):
			re = 0.0
			im = 0.0
			reWin = 0.0
			imWin = 0.0

			# for (int n = 0; n < self.__numTotalSample; n++)
			for n in range(numTotalSample):
				re = re + self.__impulseResp[n] * math.cos(2.0 * math.pi * self.__freqVec[fIndex] * n * tiSample)
				im = im - self.__impulseResp[n] * math.sin(2.0 * math.pi * self.__freqVec[fIndex] * n * tiSample)
				reWin = reWin + self.__windowedImpResp[n] * math.cos(2.0 * math.pi * self.__freqVec[fIndex] * n * tiSample)
				imWin = imWin - self.__windowedImpResp[n] * math.sin(2.0 * math.pi * self.__freqVec[fIndex] * n * tiSample)

			self.__impRespMag[fIndex] = 10.0 * math.log10(re * re + im * im)
			self.__winRespMag[fIndex] = 10.0 * math.log10(reWin * reWin + imWin * imWin)

	def __GetGainAtCutOff(self) -> float:
		re = 0.0
		im = 0.0

		numTotalSample = self.__numTotalSample
		# for (int n = 0; n < self.__numTotalSample; n++):
		for n in range(numTotalSample):
			re = re + self.__impulseResp[n] * math.cos(2.0 * math.pi * frqCutOffLo * n * tiSample)
			im = im - self.__impulseResp[n] * math.sin(2.0 * math.pi * frqCutOffLo * n * tiSample)

		return (10.0 * math.log10(re * re + im * im))

	def __ComputeWindowDFT(self):
		numFreqSample = self.__numFreqSample
		# winMag = new float[self.__numFreqSample]
		self.__winMag = [0.0] * numFreqSample

		numTotalSample = self.__numTotalSample
		tiSample = self.__tiSample
		# for (int fIndex = 0; fIndex < NUM_FREQ_SAMPLES; fIndex++):
		for fIndex in range(numFreqSample):

			re = 0.0
			im = 0.0

			# for (int n = 0; n < NUM_TOTAL_SAMPLES; n++):
			for n in range(numTotalSample):
				re = re + self.__window[n] * math.cos(2.0 * math.pi * self.__freqVec[fIndex] * n * tiSample)
				im = im - self.__window[n] * math.sin(2.0 * math.pi * self.__freqVec[fIndex] * n * tiSample)

			self.__winMag[fIndex] = 10.0 * math.log10(re * re + im * im)

	def __UpdatePlotSettings(self):

		stRespPltTyp = RespPltType(self.__stRespPltTypVar.get())
		pv(stRespPltTyp)
		stDisplay = self.__stDisplayVar.get()
		pv(stDisplay)

		# xMin1 = min(self.__dataFilrTiDomain.lines[0].get_xdata())
		# xMax1 = max(self.__dataFilrTiDomain.lines[0].get_xdata())
		yMin1 = yMin2 = yMax1 = yMax2 = 0.0

		if (stRespPltTyp == RespPltType.Impulse):		# Impulse
			if (stDisplay == "Both"):
				self.__dataFilrTiDomain.lines[0].set_visible(True)
				self.__dataFilrTiDomain.lines[1].set_visible(True)
				yMin1 = min(min(self.__dataFilrTiDomain.lines[0].get_ydata()), min(self.__dataFilrTiDomain.lines[1].get_ydata()))
				yMax1 = max(max(self.__dataFilrTiDomain.lines[0].get_ydata()), max(self.__dataFilrTiDomain.lines[1].get_ydata()))
			elif (stDisplay == "Raw Only"):
				self.__dataFilrTiDomain.lines[0].set_visible(True)
				self.__dataFilrTiDomain.lines[1].set_visible(False)
				yMin1 = min(self.__dataFilrTiDomain.lines[0].get_ydata())
				yMax1 = max(self.__dataFilrTiDomain.lines[0].get_ydata())
			elif (stDisplay == "Win Only"):
				self.__dataFilrTiDomain.lines[0].set_visible(False)
				self.__dataFilrTiDomain.lines[1].set_visible(True)
				yMin1 = min(self.__dataFilrTiDomain.lines[1].get_ydata())
				yMax1 = max(self.__dataFilrTiDomain.lines[1].get_ydata())
			self.__dataFilrTiDomain.lines[2].set_visible(False)
			self.__dataFilrTiDomain.lines[3].set_visible(False)

		else:		# Step
			if (stDisplay == "Both"):
				self.__dataFilrTiDomain.lines[2].set_visible(True)
				self.__dataFilrTiDomain.lines[3].set_visible(True)
				yMin1 = min(min(self.__dataFilrTiDomain.lines[2].get_ydata()), min(self.__dataFilrTiDomain.lines[3].get_ydata()))
				yMax1 = max(max(self.__dataFilrTiDomain.lines[2].get_ydata()), max(self.__dataFilrTiDomain.lines[3].get_ydata()))
			elif (stDisplay == "Raw Only"):
				self.__dataFilrTiDomain.lines[2].set_visible(True)
				self.__dataFilrTiDomain.lines[3].set_visible(False)
				yMin1 = min(self.__dataFilrTiDomain.lines[2].get_ydata())
				yMax1 = max(self.__dataFilrTiDomain.lines[2].get_ydata())
			elif (stDisplay == "Win Only"):
				self.__dataFilrTiDomain.lines[2].set_visible(False)
				self.__dataFilrTiDomain.lines[3].set_visible(True)
				yMin1 = min(self.__dataFilrTiDomain.lines[3].get_ydata())
				yMax1 = max(self.__dataFilrTiDomain.lines[3].get_ydata())
			self.__dataFilrTiDomain.lines[0].set_visible(False)
			self.__dataFilrTiDomain.lines[1].set_visible(False)

		if (stDisplay == "Both"):
			self.__dataFilrFreqDomain.lines[0].set_visible(True)
			self.__dataFilrFreqDomain.lines[1].set_visible(True)
			yMin2 = min(min(self.__dataFilrFreqDomain.lines[0].get_ydata()), min(self.__dataFilrFreqDomain.lines[1].get_ydata()))
			yMax2 = max(max(self.__dataFilrFreqDomain.lines[0].get_ydata()), max(self.__dataFilrFreqDomain.lines[1].get_ydata()))
		elif (stDisplay == "Raw Only"):
			self.__dataFilrFreqDomain.lines[0].set_visible(True)
			self.__dataFilrFreqDomain.lines[1].set_visible(False)
			yMin2 = min(self.__dataFilrFreqDomain.lines[0].get_ydata())
			yMax2 = max(self.__dataFilrFreqDomain.lines[0].get_ydata())
		elif (stDisplay == "Win Only"):
			self.__dataFilrFreqDomain.lines[0].set_visible(False)
			self.__dataFilrFreqDomain.lines[1].set_visible(True)
			yMin2 = min(self.__dataFilrFreqDomain.lines[1].get_ydata())
			yMax2 = max(self.__dataFilrFreqDomain.lines[1].get_ydata())

		# chrtFilrTiDomain.ChartAreas[0].RecalculateAxesScale()
		# self.__dataFilrTiDomain.ax.set_xlim(xMin1, xMax1)
		self.__dataFilrTiDomain.ax.set_ylim(yMin1, yMax1)
		self.__dataFilrTiDomain.canvas.draw()

		# self.__dataFilrFreqDomain.ax.set_xlim(xMin2, xMax2)
		self.__dataFilrFreqDomain.ax.set_ylim(yMin2, yMax2)
		self.__dataFilrFreqDomain.canvas.draw()

	def __radViewImpulse_CheckedChanged(self):
		self.__UpdatePlotSettings()

	def __cmbDisplay_SelectedIndexChanged(self, event):
		self.__UpdatePlotSettings()

	def __btnDesignFilr_Click(self):

		self.__tiSample = self.__tiSampleVar.get()
		pv(self.__tiSample)
		self.__numTotalSample = self.__numTotalSampleVar.get()
		pv(self.__numTotalSample)
		self.__numShftSample = self.__numShftSampleVar.get()
		pv(self.__numShftSample)
		self.__stWinTyp = self.__stWinTypVar.get()
		pv(self.__stWinTyp)
		self.__stFiltTyp = FilrType(self.__stFiltTypVar.get())
		pv(self.__stFiltTyp)

		if (self.__stFiltTyp == FilrType.LowPass):
			self.__frqCutOffLo = self.__frqCutOffLoVar.get()
			pv(self.__frqCutOffLo)
		elif (self.__stFiltTyp == FilrType.HighPass):
			self.__frqCutOffHi = self.__frqCutOffHiVar.get()
			pv(self.__frqCutOffHi)
		elif (self.__stFiltTyp == FilrType.BandPass):
			self.__frqCutOffLo = self.__frqCutOffLoVar.get()
			self.__frqCutOffHi = self.__frqCutOffHiVar.get()
			pv(self.__frqCutOffLo)
			pv(self.__frqCutOffHi)
		elif (self.__stFiltTyp == FilrType.BandStop):
			self.__frqCutOffLo = self.__frqCutOffLoVar.get()
			self.__frqCutOffHi = self.__frqCutOffHiVar.get()
			pv(self.__frqCutOffLo)
			pv(self.__frqCutOffHi)

		if (self.__tiSample < 0.0):
			self.__ShowErr("Sampling frequency cannot be negative.")
			return

		if (self.__frqCutOffLo >= 0.5 / self.__tiSample or self.__frqCutOffHi >= 0.5 / self.__tiSample):
			self.__ShowErr("Cut-off frequency has to be less than the Nyquist frequency (i.e. sampling freq / 2).")
			return

		if (self.__numTotalSample < 0 or self.__numShftSample < 0):
			self.__ShowErr("Total number of samples and sample shift number both need to be integers, greater than zero.")
			return

		self.__DesignFilter()

	def __cmbWin_SelectedIndexChanged(self, event):

		# if (cmbWin.SelectedIndex >= 0):

		# self.__stWinTyp = (WinType)Convert.ToInt32(cmbWin.SelectedIndex)
		self.__stWinTyp = self.__stWinTypVar.get()
		pv(self.__stWinTyp)

		# self.__ComputeTimeVec()
		self.__ComputeFreqVec()
		self.__ComputeWindow()
		self.__ComputeWindowDFT()

		'''
		chrtWinTiDomain.lines[0].Points.DataBindXY(self.__timeVec, self.__window)

		chrtWinFreqDomain.lines[0].Points.DataBindXY(self.__freqVec, self.__winMag)

		chrtWinTiDomain.ChartAreas[0].RecalculateAxesScale()
		chrtWinFreqDomain.ChartAreas[0].RecalculateAxesScale()
		'''
		self.__dataWinTiDomain.lines[0].set_data(self.__timeVec, self.__window)
		self.__dataWinTiDomain.ax.set_xlim(self.__timeVec[0], self.__timeVec[-1])
		self.__dataWinTiDomain.ax.set_ylim(min(self.__window), max(self.__window))
		self.__dataWinTiDomain.canvas.draw()
		self.__dataWinFreqDomain.lines[0].set_data(self.__freqVec, self.__winMag)
		self.__dataWinFreqDomain.ax.set_xlim(self.__freqVec[0], self.__freqVec[-1])
		self.__dataWinFreqDomain.ax.set_ylim(min(self.__winMag), max(self.__winMag))
		self.__dataWinFreqDomain.canvas.draw()

	def __infoToolStripMenuItem_Click(self):

		# MessageBox.Show("FIR Filter Designer\nWritten by Philip M. Salmony\n29 November 2019\nphilsal.co.uk", "Info", MessageBoxButtons.OK, MessageBoxIcon.Information)
		self.__ShowInfo("Info", "FIR Filter Designer\nWritten by Philip M. Salmony\n29 November 2019\nphilsal.co.uk")

	def __FilrTypeChange(self):
		self.__stFiltTyp = FilrType(self.__stFiltTypVar.get())
		pv(self.__stFiltTyp)

		if (self.__stFiltTyp == FilrType.LowPass):
			self.__txtCutOffFrqLo.configure(state='normal')
			self.__txtCutOffFrqHi.configure(state='disabled')
		elif (self.__stFiltTyp == FilrType.HighPass):
			self.__txtCutOffFrqLo.configure(state='disabled')
			self.__txtCutOffFrqHi.configure(state='normal')
		elif (self.__stFiltTyp == FilrType.BandPass):
			self.__txtCutOffFrqLo.configure(state='normal')
			self.__txtCutOffFrqHi.configure(state='normal')
		elif (self.__stFiltTyp == FilrType.BandStop):
			self.__txtCutOffFrqLo.configure(state='normal')
			self.__txtCutOffFrqHi.configure(state='normal')

	def __radPS_CheckedChanged(self):
		self.__FilrTypeChange()

	def __exportCoefficientsToolStripMenuItem_Click(self):

		filetypes = [('Text File', '*.txt'), ('All File', '*.*')]
		title  = "Export Filter Coefficients"

		r = tkFileDialog.asksaveasfilename(title=title, filetypes=filetypes)
		pv(r)

		if (r != ""):

			ext = os.path.splitext(r)[-1]
			if not ext:
				r += ".txt"

			numTotalSample = self.__numTotalSample
			with open(r, 'w') as file:
				# string[] data = new string[3]
				dataLst = [""] * 3
				# dataLst[0] = "Filr Order: " + self.__numTotalSample + " Sampling Freq (Hz): " + (1.0 / self.__tiSample).ToString("F6") + " Cut-Off Freq Lo (Hz): " + self.__frqCutOffLo.ToString("F6") + " Cut-Off Freq Hi (Hz): " + self.__frqCutOffHi.ToString("F6") + "\n\n"
				dataLst[0] = f"Filter Order: {numTotalSample} Sampling Frequency: {1.0 / self.__tiSample}Hz, Cut-Off Frequency Low: {self.__frqCutOffLo}Hz, Cut-Off Frequency High: {self.__frqCutOffHi}Hz\n\n"

				dataLst[1] = f"Coefficient, {self.__windowedImpResp[0]}"	# F7
				dataLst[2] = "float Coefficient[] = {" + str(self.__windowedImpResp[0]) + "f"	# F7

				# for (int n = 1 n < NUM_TOTAL_SAMPLES n++)
				for n in range(numTotalSample):
					dataLst[1] += "," + str(self.__windowedImpResp[n])			# F9
					dataLst[2] += "," + str(self.__windowedImpResp[n]) + "f"	# F7

				# dataLst[1] = ", ".join()
                # dataLst[2]

				dataLst[1] += "\n\n"
				dataLst[2] += "}"

				# System.IO.File.WriteAllLines(saveFileDialog.FileName, data)
				for data in dataLst:
					file.write(data)

				# MessageBox.Show("Coefficients written to file!", "Export Coefficients", MessageBoxButtons.OK, MessageBoxIcon.Information)
				self.__ShowInfo(title, f"Coefficients written to {r}")


	def __exportTiDomainDataToolStripMenuItem_Click(self):

		filetypes = [('Text File', '*.txt'), ('All File', '*.*')]
		title  = "Export Time Domain Data"

		r = tkFileDialog.asksaveasfilename(title=title, filetypes=filetypes)
		pv(r)

		if (r != ""):
			ext = os.path.splitext(r)[-1]
			if not ext:
				r += ".txt"

			numTotalSample = self.__numTotalSample
			with open(r, 'w') as file:
				# string[] data = new string[4]
				dataLst = [""] * 4
				# dataLst[0] = "[TIME DOMAIN DATA (TIME/IMPULSE/STEP)] Filter Order: " + self.__numTotalSample + " Sampling Frequency (Hz): " + (1.0 / self.__tiSample).ToString("F6") + " Cut-Off Frequency Low (Hz): " + self.__frqCutOffLo.ToString("F6") + " Cut-Off Frequency High (Hz): " + self.__frqCutOffHi.ToString("F6") + "\n\n"
				dataLst[0] = f"[TIME DOMAIN DATA (TIME/IMPULSE/STEP)]\n\nFilter Order: {numTotalSample}, Sampling Frequency: {1.0 / self.__tiSample}Hz, Cut-Off Frequency Low: {self.__frqCutOffLo}Hz, Cut-Off Frequency High: {self.__frqCutOffHi}Hz\n\n"

				dataLst[1] = f"Time Vector (s), {self.__timeVec[0]}"				# F6
				dataLst[2] = f"Windowed Impulse Response, {self.__windowedImpResp[0]}"  # F9
				dataLst[3] = f"Windowed Step Response, {self.__windowedStepResp[0]}"     # F9
				# for (int n = 1; n < NUM_TOTAL_SAMPLES; n++):
				for n in range(numTotalSample):
					dataLst[1] += "," + str(self.__timeVec[n])              	# F6
					dataLst[2] += "," + str(self.__windowedImpResp[n])    # F9
					dataLst[3] += "," + str(self.__windowedStepResp[n])       # F9
				# dataLst[1] = ", ".join(str(self.__timeVec))
				# dataLst[2] = ", ".join(str(self.__windowedImpResp))
				# dataLst[3] = ", ".join(str(self.__windowedStepResp))

				dataLst[1] += "\n\n"
				dataLst[2] += "\n\n"

				# System.IO.File.WriteAllLines(saveFileDialog.FileName, data)
				for data in dataLst:
					file.write(data)
				# MessageBox.Show("Data written to file!", "Export Time Domain Data", MessageBoxButtons.OK, MessageBoxIcon.Information)
				self.__ShowInfo(title, f"Data written to {r}")


	def __exportFreqDomainDataToolStripMenuItem_Click(self):

		filetypes = [('Text File', '*.txt'), ('All File', '*.*')]
		title  = "Export Frequency Domain Data"

		r = tkFileDialog.asksaveasfilename(title=title, filetypes=filetypes)
		pv(r)

		if (r != ""):
			ext = os.path.splitext(r)[-1]
			if not ext:
				r += ".txt"

			with open(r, 'w') as file:

				# string[] data = new string[4]
				dataLst = [""] * 4
				# dataLst[0] = "[FREQUENCY DOMAIN DATA (FREQUENCY/RAW/WINDOWED)] Filr Order: " + self.__numTotalSample + " Sampling Freq (Hz): " + (1.0 / self.__tiSample).ToString("F6") + " Cut-Off Freq Lo (Hz): " + self.__frqCutOffLo.ToString("F6") + " Cut-Off Freq Hi (Hz): " + self.__frqCutOffHi.ToString("F6") + "\n\n"
				dataLst[0] = f"[FREQUENCY DOMAIN DATA (FREQUENCY/RAW/WINDOWED)]\n\nFilter Order: {self.__numTotalSample}, Sampling Frequency: {1.0 / self.__tiSample}Hz, Cut-Off Frequency Low: {self.__frqCutOffLo}Hz, Cut-Off Frequency High: {self.__frqCutOffHi}Hz\n\n"

				dataLst[1] = f"Frequency Vecctor (Hz), {self.__freqVec[0]}"		# F6
				dataLst[2] = f"Impulse Response Magnitude, {self.__impRespMag[0]}"		# F9
				dataLst[3] = f"Window Response Magnitude, {self.__winRespMag[0]}" 	# F9
				# for (int n = 1; n < NUM_FREQ_SAMPLES; n++):
				for n in range(self.__numFreqSample):
					dataLst[1] += "," + str(self.__freqVec[n])			# F6
					dataLst[2] += "," + str(self.__impRespMag[n])    # F9
					dataLst[3] += "," + str(self.__winRespMag[n])    # F9

				dataLst[1] += "\n\n"
				dataLst[2] += "\n\n"

				# System.IO.File.WriteAllLines(saveFileDialog.FileName, data)
				for data in dataLst:
					file.write(data)
				# MessageBox.Show("Data written to file!", "Export Frequency Domain Data", MessageBoxButtons.OK, MessageBoxIcon.Information)
				self.__ShowInfo(title, f"Data written to {r}")

if __name__ == '__main__':
	root = tk.Tk()
	myapp = App(root)
	myapp.mainloop()