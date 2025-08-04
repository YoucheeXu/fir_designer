#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import math
from enum import IntEnum, unique, Enum

# from tkinter import *
import tkinter as tk
import tkinter.messagebox as tkMessageBox
import tkinter.filedialog as tkFileDialog
import tkinter.ttk as ttk

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

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
class FilterType(Enum):
	LowPass = 0
	HighPass = 1
	BandPass = 2
	BandStop = 3

class MatPlotData:
    def __init__(self):
        self.canvas = None
        self.ax = None

class App(tk.Frame):

	def __init__(self, master):
		super().__init__(master)

		# Initial settings
		
		# Parameters
		self.__SAMPLE_TIME_S = 0.01
		self.__CUTOFF_FREQLO_HZ = 20.0
		self.__CUTOFF_FREQHI_HZ = 20.0
		self.__NUM_TOTAL_SAMPLES = 64
		self.__NUM_SHIFT_SAMPLES = 32
		self.__WIN_TYPE = WinType.Hamming
		self.__FILT_TYPE = FilterType.LowPass
		self.__NUM_FREQ_SAMPLES = 256

		# Time Domain
		self.__timeVec = []
		self.__impulseResp = []
		self.__stepResp = []
		self.__window = []
		self.__windowedImpulseResp = []
		self.__windowedStepResp = []

		# Freq Domain
		self.__freqVecHz = []
		self.__impRespMag = []
		self.__winRespMag = []
		self.__winMag = []

		# cmbWin.SelectedIndex = 5
		# radViewImpulse.Checked = True
		# cmbDisplay.SelectedIndex = 0

		self.__CreateWindow(master)

		self.__DesignFilter()

	def __ShowErr(self, err):
		tkMessageBox.showerror("Design FIR", err)

	def __CenterWindow(self, window, width, hight):
		'''
			设置窗口居中和宽高
			:param window:主窗体
			:param Width:窗口宽度
			:param Hight:窗口高度
			:return:无
		'''
		# 获取屏幕宽度和高度
		sw = window.winfo_screenwidth()
		sh = window.winfo_screenheight()

		# 计算中心坐标
		cen_x = (sw - width) / 2
		cen_y = (sh - hight) / 2

		# 设置窗口大小并居中
		window.geometry('%dx%d+%d+%d' % (width, hight, cen_x, cen_y))

	def __CreateWindow(self, root):

		root.title("Design FIR")  # 设置窗口标题
		self.__CenterWindow(root, 1000, 800)

		framTimeDomain = tk.Frame(root)

		# chrtFilterTimeDomain = tk.LabelFrame(framTimeDomain, text="Filter Time Domain")
		self.__dataFilterTimeDomain = MatPlotData()
		fig, self.__dataFilterTimeDomain.ax = self.__CreateMatplot(r"Filter Time Domain", r"Time [s]", r"Time Domain")
		self.__dataFilterTimeDomain.canvas = FigureCanvasTkAgg(fig, master=framTimeDomain)  # A tk.DrawingArea.
		self.__dataFilterTimeDomain.canvas.draw()
		self.__dataFilterTimeDomain.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
		# chrtFilterTimeDomain.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

		fig, ax = self.__CreateMatplot(r"Window Time Domain", r"Time [s]", r"Window Time Domain")
		canvasWinTimeDomain = FigureCanvasTkAgg(fig, master=framTimeDomain)  # A tk.DrawingArea.
		canvasWinTimeDomain.draw()
		canvasWinTimeDomain.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

		framTimeDomain.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

		framFreqDomain = tk.Frame(root)

		fig, ax = self.__CreateMatplot(r"Filter Frequency Domain", r"Frequency [Hz]", r"Filter Frequency Domain")
		canvasFilterFreqDomain = FigureCanvasTkAgg(fig, master=framFreqDomain)  # A tk.DrawingArea.
		canvasFilterFreqDomain.draw()
		canvasFilterFreqDomain.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

		fig, ax = self.__CreateMatplot(r"Window Frequency Domain", r"Frequency [Hz]", r"Window Frequency Domain")
		canvasWinFreqDomain = FigureCanvasTkAgg(fig, master=framFreqDomain)  # A tk.DrawingArea.
		canvasWinFreqDomain.draw()
		canvasWinFreqDomain.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

		framFreqDomain.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

		frmParas = tk.LabelFrame(root, text="Parameters")

		# radBS radBP radHP radLP
		self.__radFiltTyp = tk.IntVar()
		self.__radFiltTyp.set(0)	# 初始化变量值
		frmFiltTyp = tk.LabelFrame(frmParas, text="Filter Type")
		radio = ttk.Radiobutton(frmFiltTyp, variable = self.__radFiltTyp, value = 0, text ='Low Pass', command = self.__radPS_CheckedChanged)
		radio.grid(row=0, column=0)
		radio = ttk.Radiobutton(frmFiltTyp, variable = self.__radFiltTyp, value = 1, text ='High Pass', command = self.__radPS_CheckedChanged)
		radio.grid(row=0, column=1)
		radio = ttk.Radiobutton(frmFiltTyp, variable = self.__radFiltTyp, value = 2, text ='Band Pass', command = self.__radPS_CheckedChanged)
		radio.grid(row=1, column=0)
		radio = ttk.Radiobutton(frmFiltTyp, variable = self.__radFiltTyp, value = 3, text ='Band Stop', command = self.__radPS_CheckedChanged)
		radio.grid(row=1, column=1)
		frmFiltTyp.grid(row=0, column=4)

		# cmbDisplay
		
		# radViewStep radViewImpulse

		self.__cmbWinTyp = tk.StringVar()
		self.__cmbWinTyp.set("0")
		lblcmbWinTyp = tk.Label(frmParas, text = "Type of Win: ")  
		lblcmbWinTyp.grid(row=2, column=2)
		options = ["Rectangular", "Triangular", "Welch", "Sine", "Hann", "Hamming", "Blackman", "Nuttall", "BlackmanNuttall", "BlackmanHarris", "FlatTop"]
		cmbcmbWinTyp = ttk.Combobox(frmParas, textvariable = self.__cmbWinTyp, values = options, state='readonly', width=10)
		cmbcmbWinTyp.current(int(self.__WIN_TYPE))
		cmbcmbWinTyp.bind('<<ComboboxSelected>>', self.__cmbWin_SelectedIndexChanged)
		cmbcmbWinTyp.grid(row=2, column=3)

		self._txtCufOffFreqLo = tk.DoubleVar()
		self._txtCufOffFreqLo.set(self.__CUTOFF_FREQLO_HZ)
		lblCutOffFreq = tk.Label(frmParas, text = "Low Frquency(Hz):")  
		lblCutOffFreq.grid(row=0, column=2)
		enyCutOffFreq = ttk.Entry(frmParas, textvariable = self._txtCufOffFreqLo, width=10)
		enyCutOffFreq.grid(row=0, column=3)

		self._txtCufOffFreqHi = tk.DoubleVar()
		self._txtCufOffFreqHi.set(self.__CUTOFF_FREQHI_HZ)
		lblCufOffFrqHi = tk.Label(frmParas, text = "High Frquency(Hz):")  
		lblCufOffFrqHi.grid(row=1, column=2)
		enyCufOffFrqHi = ttk.Entry(frmParas, textvariable = self._txtCufOffFreqHi, width=10)
		enyCufOffFrqHi.grid(row=1, column=3)

		self.__txtShiftSamples = tk.IntVar()
		self.__txtShiftSamples.set(self.__NUM_SHIFT_SAMPLES)
		lblShiftSamples = tk.Label(frmParas, text = "Shift Samples:")  
		lblShiftSamples.grid(row=0, column=0)
		enyShiftSamples = ttk.Entry(frmParas, textvariable = self.__txtShiftSamples, width=10)
		enyShiftSamples.grid(row=0, column=1)

		self.__txtFilterLength = tk.IntVar()
		self.__txtFilterLength.set(self.__NUM_TOTAL_SAMPLES)
		lblFilterLength = tk.Label(frmParas, text = "Length of Filter:")  
		lblFilterLength.grid(row=1, column=0)
		enyFilterLength = ttk.Entry(frmParas, textvariable = self.__txtFilterLength, width=10)
		enyFilterLength.grid(row=1, column=1)

		self.__txtSampleTime = tk.DoubleVar()
		self.__txtSampleTime.set(self.__SAMPLE_TIME_S)
		lblSampFrq = tk.Label(frmParas, text = "Sample Time(s):")
		lblSampFrq.grid(row=2, column=0)
		enySampFrq = ttk.Entry(frmParas, textvariable = self.__txtSampleTime, width=10)
		enySampFrq.grid(row=2, column=1)

		frmParas.pack(side=tk.LEFT)

		btnDesignFilter = ttk.Button(root, text ="Design Filter", command = self.__btnDesignFilter_Click)
		btnDesignFilter.pack(side=tk.LEFT)

		# 创建一个顶级菜单
		menubar = tk.Menu(root)

		# 在顶层窗口menubar中加入下拉窗口filemenu
		filemenu = tk.Menu(menubar, tearoff=False)

		filemenu.add_command(label="Export Coefficients...", command=self.__exportCoefficientsToolStripMenuItem_Click)
		filemenu.add_command(label="Export Time Domain Data...", command=self.__exportTimeDomainDataToolStripMenuItem_Click)
		filemenu.add_command(label="Export Frequency Domain Data...", command=self.__exportFreqDomainDataToolStripMenuItem_Click)
		filemenu.add_separator()

		filemenu.add_command(label="Exit", command=root.destroy)
		menubar.add_cascade(label="File", menu=filemenu)	# 为filemenu命名为‘文件’

		# infoToolStripMenuItem
		menubar.add_command(label = "Info")

		# 显示菜单
		root.config(menu=menubar)

	def __CreateMatplot(self, title, xLabel, yLabel):

		fig = plt.Figure(figsize=(0, 0), dpi=100)
		t = np.arange(0, 3, .01)
		ax = fig.add_subplot()
		line, = ax.plot(t, 2 * np.sin(2 * np.pi * t))

		ax.set_title(title)
		ax.set_xlabel(xLabel)
		ax.set_ylabel(yLabel)

		return fig, ax

	def __DesignFilter(self):
		self.__ComputeTimeVec()
		self.__ComputeWin()
		self.__ComputeResps()
		self.__ComputeWinedResps()

		self.__ComputeFreqVec()
		self.__ComputeRespBode()
		self.__ComputeWinDFT()	

		self.__UpdateCharts()

	def __UpdateCharts(self):

		'''
		chrtFilterTimeDomain.Series[0].Points.DataBindXY(self.__timeVec, self.__impulseResp)
		chrtFilterTimeDomain.Series[1].Points.DataBindXY(self.__timeVec, self.__windowedImpulseResp)
		chrtFilterTimeDomain.Series[2].Points.DataBindXY(self.__timeVec, self.__stepResp)
		chrtFilterTimeDomain.Series[3].Points.DataBindXY(self.__timeVec, self.__windowedStepResp)
		chrtFilterTimeDomain.Series[2].Enabled = False
		chrtFilterTimeDomain.Series[3].Enabled = False

		chrtFilterTimeDomain.ChartAreas[0].RecalculateAxesScale()
		chrtFilterTimeDomain.Update()
		'''
		ax = self.__dataFilterTimeDomain.ax
		title = ax.get_title()
		xLabel = ax.get_xlabel()
		yLabel = ax.get_ylabel()
		ax.cla()
		ax.set_title(title)
		ax.set_xlabel(xLabel)
		ax.set_ylabel(yLabel)
		ax.plot(self.__timeVec, self.__impulseResp)
		ax.plot(self.__timeVec, self.__windowedImpulseResp)
		# self.__dataFilterTimeDomain.ax.plot(self.__timeVec, self.__stepResp)
		# self.__dataFilterTimeDomain.ax.plot(self.__timeVec, self.__windowedStepResp)
		self.__dataFilterTimeDomain.canvas.draw()

		'''
		chrtWinTimeDomain.Series[0].Points.DataBindXY(self.__timeVec, self.__window)
		chrtWinTimeDomain.ChartAreas[0].RecalculateAxesScale()
		'''

		'''
		chrtFilterFreqDomain.Series[0].Points.DataBindXY(self.__freqVecHz, self.__impRespMag)
		chrtFilterFreqDomain.Series[1].Points.DataBindXY(self.__freqVecHz, self.__winRespMag)
		chrtFilterFreqDomain.ChartAreas[0].AxisX.Interval = 10.0 * (0.5 / self.__SAMPLE_TIME_S) / (self.__NUM_FREQ_SAMPLES - 1.0)
		chrtFilterFreqDomain.ChartAreas[0].RecalculateAxesScale()
		'''

		'''
		chrtWinFreqDomain.Series[0].Points.DataBindXY(self.__freqVecHz, self.__winMag)
		chrtWinFreqDomain.ChartAreas[0].AxisX.Interval = 10.0 * (0.5 / self.__SAMPLE_TIME_S) / (self.__NUM_FREQ_SAMPLES - 1.0)
		chrtWinFreqDomain.ChartAreas[0].RecalculateAxesScale()
		'''


	# Time Domain Functions
	def __ComputeTimeVec(self):
		# self.__timeVec = new float[self.__NUM_TOTAL_SAMPLES]
		self.__timeVec.clear()

		# for (int n = 0 n < self.__NUM_TOTAL_SAMPLES n++):
		for n in range(self.__NUM_TOTAL_SAMPLES):
			self.__timeVec.append(n * self.__SAMPLE_TIME_S)

	def __ComputeResps(self):

		# self.__impulseResp = new float[self.__NUM_TOTAL_SAMPLES]
		self.__impulseResp.clear()
		# self.__stepResp = new float[self.__NUM_TOTAL_SAMPLES]
		self.__stepResp.clear()

		# for (int n = 0 n < self.__NUM_TOTAL_SAMPLES n++)
		for n in range(self.__NUM_TOTAL_SAMPLES):
			if (n != self.__NUM_SHIFT_SAMPLES):
				if self.__FILT_TYPE == FilterType.LowPass:
					self.__impulseResp.append(math.sin(2.0 * math.pi * self.__CUTOFF_FREQLO_HZ * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)) / (math.pi * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)))
				elif self.__FILT_TYPE == FilterType.HighPass:                            
					self.__impulseResp.append((math.sin(math.pi * (n - self.__NUM_SHIFT_SAMPLES)) - math.sin(2.0 * math.pi * self.__CUTOFF_FREQLO_HZ * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES))) / (math.pi * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)))
				elif self.__FILT_TYPE == FilterType.BandPass:
					self.__impulseResp.append((math.sin(2.0 * math.pi * self.__CUTOFF_FREQHI_HZ * self.__SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES)) - math.sin(2.0 * math.pi * self.__CUTOFF_FREQLO_HZ * self.__SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES))) / (math.pi * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)))
				elif self.__FILT_TYPE == FilterType.BandStop:
					self.__impulseResp.append((math.sin(2.0 * math.pi * self.__CUTOFF_FREQLO_HZ * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)) - math.sin(2.0 * math.pi * self.__CUTOFF_FREQHI_HZ * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)) + math.sin(math.pi * (n - self.__NUM_SHIFT_SAMPLES))) / (math.pi * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)))

			else: # Avoid divide-by-zero, limit is 2*fc
				if self.__FILT_TYPE == FilterType.LowPass:
					self.__impulseResp.append(2.0 * self.__CUTOFF_FREQLO_HZ)
				elif self.__FILT_TYPE == FilterType.HighPass:
					self.__impulseResp.append(1.0 / self.__SAMPLE_TIME_S - 2.0 * self.__CUTOFF_FREQLO_HZ)
				elif self.__FILT_TYPE == FilterType.BandPass:
					self.__impulseResp.append(2.0 * self.__CUTOFF_FREQHI_HZ - 2.0 * self.__CUTOFF_FREQLO_HZ)
				elif self.__FILT_TYPE == FilterType.BandStop:
					self.__impulseResp.append(2.0 * self.__CUTOFF_FREQLO_HZ - 2.0 * self.__CUTOFF_FREQHI_HZ + 1.0 / self.__SAMPLE_TIME_S)

		# Normalise by DC gain to achieve 0dB gain at DC and then compute step response
		# for (int n = 0 n < self.__NUM_TOTAL_SAMPLES n++)
		for n in range(self.__NUM_TOTAL_SAMPLES):
			self.__impulseResp[n] *= self.__SAMPLE_TIME_S

			if (n == 0):
				self.__stepResp.append(0.5 * self.__impulseResp[n])
			else:		
				self.__stepResp.append(self.__stepResp[n - 1] + 0.5 * (self.__impulseResp[n] + self.__impulseResp[n - 1]))

	def __ComputeWin(self):

		# self.__window = new float[self.__NUM_TOTAL_SAMPLES]
		self.__window.clear()

		# for (int n = 0 n < self.__NUM_TOTAL_SAMPLES n++)
		for n in range(self.__NUM_TOTAL_SAMPLES):
			if self.__WIN_TYPE == "Rectangular":
				self.__window.append(1.0)
			elif self.__WIN_TYPE == "Triangular":
				self.__window.append(1.0 - math.abs((n - 0.5 * self.__NUM_TOTAL_SAMPLES) / (0.5 * self.__NUM_TOTAL_SAMPLES)))
			elif self.__WIN_TYPE == "Welch":
				self.__window.append(1.0 - math.pow((n - 0.5 * self.__NUM_TOTAL_SAMPLES) / (0.5 * self.__NUM_TOTAL_SAMPLES), 2.0))
			elif self.__WIN_TYPE == "Sine":
				self.__window.append(math.sin(math.pi * n / (self.__NUM_TOTAL_SAMPLES)))
			elif self.__WIN_TYPE == "Hann":
				self.__window.append(0.5 * (1 - math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES))))
			elif self.__WIN_TYPE == "Hamming":
				self.__window.append((25.0 / 46.0) - (21.0 / 46.0) * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)))
			elif self.__WIN_TYPE == "Blackman":
				self.__window.append(0.42 - 0.5 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.08 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)))
			elif self.__WIN_TYPE == "Nuttall":
				self.__window.append(0.355768 - 0.487396 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.144232 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) - 0.012604 * math.cos(6.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)))
			elif self.__WIN_TYPE == "BlackmanNuttall":
				self.__window.append(0.3635819 - 0.4891775 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.1365995 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) - 0.0106411 * math.cos(6.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)))
			elif self.__WIN_TYPE == "BlackmanHarris":
				self.__window.append(0.35875 - 0.48829 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.14128 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) - 0.01168 * math.cos(6.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)))
			elif self.__WIN_TYPE == "FlatTop":
				self.__window.append(0.21557895 - 0.41663158 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.277263158 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) - 0.083578947 * math.cos(6.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.006947368 * math.cos(8.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)))
			else:
				self.__window.append(1.0)

	def __ComputeWinedResps(self):

		# self.__windowedImpulseResp = new float[self.__NUM_TOTAL_SAMPLES]
		self.__windowedImpulseResp.clear()
		# self.__windowedStepResp = new float[self.__NUM_TOTAL_SAMPLES]
		self.__windowedStepResp = [0.0] * self.__NUM_TOTAL_SAMPLES

		# for (int n = 0; n < self.__NUM_TOTAL_SAMPLES; n++):
		for n in range(self.__NUM_TOTAL_SAMPLES):
			self.__windowedImpulseResp.append(self.__impulseResp[n] * self.__window[n])

			if (n == 0):
				self.__windowedStepResp.append(0.5 * self.__windowedStepResp[n])
			else:		
				self.__windowedStepResp.append(self.__windowedStepResp[n - 1] + 0.5 * (self.__windowedImpulseResp[n] + self.__windowedImpulseResp[n - 1]))

	# Freq Domain Functions
	def __ComputeFreqVec(self):

		# self.__freqVecHz = new float[self.__NUM_FREQ_SAMPLES]
		self.__freqVecHz.clear()

		df = (0.5 / self.__SAMPLE_TIME_S) / (self.__NUM_FREQ_SAMPLES - 1.0)

		# for (int n = 0; n < self.__NUM_FREQ_SAMPLES; n++):
		for n in range(self.__NUM_TOTAL_SAMPLES):
			self.__freqVecHz.append(n * df)

	def __ComputeRespBode(self):

		# self.__impRespMag = new float[self.__NUM_FREQ_SAMPLES]
		self.__impRespMag.clear()
		# self.__winRespMag = new float[self.__NUM_FREQ_SAMPLES]
		self.__winRespMag.clear()

		# for (int fIndex = 0; fIndex < NUM_FREQ_SAMPLES; fIndex++):
		for fIndex in range(self.__NUM_FREQ_SAMPLES):
			re = 0.0
			im = 0.0
			reWin = 0.0
			imWin = 0.0

			# for (int n = 0; n < self.__NUM_TOTAL_SAMPLES; n++)
			for n in range(self.__NUM_TOTAL_SAMPLES):
				pv(n)
				pv(self.__impulseResp[n])
				pv(self.__freqVecHz[fIndex])
				re += self.__impulseResp[n] * math.cos(2.0 * math.pi * self.__freqVecHz[fIndex] * n * self.__SAMPLE_TIME_S)
				im -= self.__impulseResp[n] * math.sin(2.0 * math.pi * self.__freqVecHz[fIndex] * n * self.__SAMPLE_TIME_S)
				reWin += self.__windowedImpulseResp[n] * math.cos(2.0 * math.pi * self.__freqVecHz[fIndex] * n * self.__SAMPLE_TIME_S)
				imWin -= self.__windowedImpulseResp[n] * math.sin(2.0 * math.pi * self.__freqVecHz[fIndex] * n * self.__SAMPLE_TIME_S)

			self.__impRespMag.append(10.0 * math.log10(re * re + im * im))
			self.__winRespMag.append(10.0 * math.log10(reWin * reWin + imWin * imWin))

	def __GetGainAtCutOff(self) -> float: 
		re = 0.0
		im = 0.0

		# for (int n = 0; n < self.__NUM_TOTAL_SAMPLES; n++):
		for n in range(self.__NUM_TOTAL_SAMPLES):
			re = re + self.__impulseResp[n] * math.cos(2.0 * math.pi * self.__CUTOFF_FREQLO_HZ * n * self.__SAMPLE_TIME_S)
			im = im - self.__impulseResp[n] * math.sin(2.0 * math.pi * self.__CUTOFF_FREQLO_HZ * n * self.__SAMPLE_TIME_S)

		return (10.0 * math.log10(re * re + im * im))

	def __ComputeWinDFT(self):

		# self.__winMag = new float[self.__NUM_FREQ_SAMPLES]
		self.__winMag.clear()

		# for (int fIndex = 0; fIndex < NUM_FREQ_SAMPLES; fIndex++):
		for fIndex in range(self.__NUM_FREQ_SAMPLES):

			re = 0.0
			im = 0.0

			# for (int n = 0; n < NUM_TOTAL_SAMPLES; n++):
			for n in range(self.__NUM_TOTAL_SAMPLES):
				re = re + self.__window[n] * math.cos(2.0 * math.pi * self.__freqVecHz[fIndex] * n * self.__SAMPLE_TIME_S)
				im = im - self.__window[n] * math.sin(2.0 * math.pi * self.__freqVecHz[fIndex] * n * self.__SAMPLE_TIME_S)

			self.__winMag.append(10.0 * math.log10(re * re + im * im))

	def __UpdatePlotSettings(self):

		if (radViewImpulse.Checked):	
			if (cmbDisplay.SelectedIndex == 0):		
				chrtFilterTimeDomain.Series[0].Enabled = True
				chrtFilterTimeDomain.Series[1].Enabled = True
			elif (cmbDisplay.SelectedIndex == 1):		
				chrtFilterTimeDomain.Series[0].Enabled = True
				chrtFilterTimeDomain.Series[1].Enabled = False		
			elif (cmbDisplay.SelectedIndex == 2):		
				chrtFilterTimeDomain.Series[0].Enabled = False
				chrtFilterTimeDomain.Series[1].Enabled = True
			chrtFilterTimeDomain.Series[2].Enabled = False
			chrtFilterTimeDomain.Series[3].Enabled = False

		else:
			if (cmbDisplay.SelectedIndex == 0):		
				chrtFilterTimeDomain.Series[2].Enabled = True
				chrtFilterTimeDomain.Series[3].Enabled = True		
			elif (cmbDisplay.SelectedIndex == 1):		
				chrtFilterTimeDomain.Series[2].Enabled = True
				chrtFilterTimeDomain.Series[3].Enabled = False		
			elif (cmbDisplay.SelectedIndex == 2):		
				chrtFilterTimeDomain.Series[2].Enabled = False
				chrtFilterTimeDomain.Series[3].Enabled = True
			chrtFilterTimeDomain.Series[0].Enabled = False
			chrtFilterTimeDomain.Series[1].Enabled = False

		if (cmbDisplay.SelectedIndex == 0):	
			chrtFilterFreqDomain.Series[0].Enabled = True
			chrtFilterFreqDomain.Series[1].Enabled = True	
		elif (cmbDisplay.SelectedIndex == 1):	
			chrtFilterFreqDomain.Series[0].Enabled = True
			chrtFilterFreqDomain.Series[1].Enabled = False	
		elif (cmbDisplay.SelectedIndex == 2):	
			chrtFilterFreqDomain.Series[0].Enabled = False
			chrtFilterFreqDomain.Series[1].Enabled = True

		chrtFilterTimeDomain.ChartAreas[0].RecalculateAxesScale()


	def __radViewImpulse_CheckedChanged(self):
		self.__UpdatePlotSettings()

	def __cmbDisplay_SelectedIndexChanged(self):
		self.__UpdatePlotSettings()

	def __btnDesignFilter_Click(self):

		# self.__SAMPLE_TIME_S = 1.0 / float(txtSamplingFreq.get())
		self.__SAMPLE_TIME_S = self.__txtSampleTime.get()
		self.__NUM_TOTAL_SAMPLES = self.__txtFilterLength.get()
		self.__NUM_SHIFT_SAMPLES = self.__txtShiftSamples.get()
		self.__FILT_TYPE = FilterType(self.__radFiltTyp.get())
		pv(self.__FILT_TYPE)
		self.__WIN_TYPE = self.__cmbWinTyp.get()
		pv(self.__WIN_TYPE)

		if (self.__FILT_TYPE == FilterType.LowPass):
			self.__CUTOFF_FREQLO_HZ = self._txtCufOffFreqHi.get()
		elif (self.__FILT_TYPE == FilterType.HighPass):
			self.__CUTOFF_FREQLO_HZ = self._txtCufOffFreqLo.get()
		elif (self.__FILT_TYPE == FilterType.BandPass):
			self.__CUTOFF_FREQLO_HZ = self._txtCufOffFreqLo.get()
			self.__CUTOFF_FREQHI_HZ = self._txtCufOffFreqHi.get()
		elif (self.__FILT_TYPE == FilterType.BandStop):
			self.__CUTOFF_FREQLO_HZ = self._txtCufOffFreqLo.get()
			self.__CUTOFF_FREQHI_HZ = self._txtCufOffFreqHi.get()

		if (self.__SAMPLE_TIME_S < 0.0):
			self.__ShowErr("Sampling freq cannot be negative.")
			return

		if (self.__CUTOFF_FREQLO_HZ >= 0.5 / self.__SAMPLE_TIME_S or self.__CUTOFF_FREQHI_HZ >= 0.5 / self.__SAMPLE_TIME_S):	
			self.__ShowErr("Cut-off freq has to be less than the Nyquist freq (i.e. sampling freq / 2).")
			return

		if (self.__NUM_TOTAL_SAMPLES < 0 or self.__NUM_SHIFT_SAMPLES < 0):	
			self.__ShowErr("Total number of samples and sample shift number both need to be integers, greater than zero.")
			return

		self.__DesignFilter()

		# self.__UpdatePlotSettings()

	def __cmbWin_SelectedIndexChanged(self, event):

		# if (cmbWin.SelectedIndex >= 0):

		# self.__WIN_TYPE = (WinType)Convert.ToInt32(cmbWin.SelectedIndex)
		self.__WIN_TYPE = self.__cmbWinTyp.get()
		pv(self.__WIN_TYPE)

		self.__ComputeTimeVec()
		self.__ComputeFreqVec()
		self.__ComputeWin()
		self.__ComputeWinDFT()

		'''
		chrtWinTimeDomain.Series[0].Points.DataBindXY(self.__timeVec, self.__window)

		chrtWinFreqDomain.Series[0].Points.DataBindXY(self.__freqVecHz, self.__winMag)

		chrtWinTimeDomain.ChartAreas[0].RecalculateAxesScale()
		chrtWinFreqDomain.ChartAreas[0].RecalculateAxesScale()
		'''
	
	def __infoToolStripMenuItem_Click(self):

		MessageBox.Show("FIR Filter Designer\nWritten by Philip M. Salmony\n29 November 2019\nphilsal.co.uk", "Info", MessageBoxButtons.OK, MessageBoxIcon.Information)


	def __FilterTypeChange(self):
		print(self.__radFiltTyp)

		if (radLP.Checked):	
			txtCutOffFreq.Enabled = False
			txtCutOffFreqHigh.Enabled = True
		elif (radHP.Checked):	
			txtCutOffFreq.Enabled = True
			txtCutOffFreqHigh.Enabled = False
		elif (radBP.Checked):	
			txtCutOffFreq.Enabled = True
			txtCutOffFreqHigh.Enabled = True
		elif (radBS.Checked):	
			txtCutOffFreq.Enabled = True
			txtCutOffFreqHigh.Enabled = True

	def __radPS_CheckedChanged(self):
		self.__FilterTypeChange()

	def __exportCoefficientsToolStripMenuItem_Click(self):

		filetypes = [('Text File', '*.txt')]
		title  = "Export Filter Coefficients"

		r = tkFileDialog.asksaveasfilename(title=title, filetypes=filetypes)
		pv(r)

		if (r != ""):

			ext = os.path.splitext(r)[-1]
			if not ext:
				r += ".txt"

			with open(r, 'w') as file:
				# string[] data = new string[3]
				data = [""] * 3
				# data[0] = "Filter Order: " + self.__NUM_TOTAL_SAMPLES + " Sampling Freq (Hz): " + (1.0 / self.__SAMPLE_TIME_S).ToString("F6") + " Cut-Off Freq Lo (Hz): " + self.__CUTOFF_FREQLO_HZ.ToString("F6") + " Cut-Off Freq Hi (Hz): " + self.__CUTOFF_FREQHI_HZ.ToString("F6") + "\n\n"
				data[0] = "Filter Order: " + str(self.__NUM_TOTAL_SAMPLES) + " Sampling Freq (Hz): " + str(1.0 / self.__SAMPLE_TIME_S) + " Cut-Off Freq Lo (Hz): " + str(self.__CUTOFF_FREQLO_HZ) + " Cut-Off Freq Hi (Hz): " + str(self.__CUTOFF_FREQHI_HZ) + "\n\n"

				data[1] = str(self.__windowedImpulseResp[0])	# F7
				data[2] = "float coeff[] = {" + str(self.__windowedImpulseResp[0]) + "f"	# F7
				# for (int n = 1 n < NUM_TOTAL_SAMPLES n++)
				for n in range(self.__NUM_TOTAL_SAMPLES):
					data[1] += "," + str(self.__windowedImpulseResp[n])			# F9
					data[2] += "," + str(self.__windowedImpulseResp[n]) + "f"	# F7
			
				data[1] += "\n\n"
				data[2] += "}"
				
				# System.IO.File.WriteAllLines(saveFileDialog.FileName, data)
				r.write(data)
			
			MessageBox.Show("Coefficients written to file!", "Export Coefficients", MessageBoxButtons.OK, MessageBoxIcon.Information)
	

	def __exportTimeDomainDataToolStripMenuItem_Click(self):

		saveFileDialog.Filter = "Text File|*.txt"
		saveFileDialog.Title  = "Export Time Domain Data"
		saveFileDialog.ShowDialog()

		if (saveFileDialog.FileName != ""):
	
			# string[] data = new string[4]
			data = [""] * 4
			data[0] = "[TIME DOMAIN DATA (TIME/IMPULSE/STEP)] Filter Order: " + self.__NUM_TOTAL_SAMPLES + " Sampling Freq (Hz): " + (1.0 / self.__SAMPLE_TIME_S).ToString("F6") + " Cut-Off Freq Lo (Hz): " + self.__CUTOFF_FREQLO_HZ.ToString("F6") + " Cut-Off Freq Hi (Hz): " + self.__CUTOFF_FREQHI_HZ.ToString("F6") + "\n\n"

			data[1] = self.__timeVec[0].ToString("F6")
			data[2] = self.__windowedImpulseResp[0].ToString("F9")
			data[3] = self.__windowedStepResp[0].ToString("F9")

			# for (int n = 1; n < NUM_TOTAL_SAMPLES; n++):
			for n in range(self.__NUM_TOTAL_SAMPLES):
				data[1] += "," + self.__timeVec[n].ToString("F6")
				data[2] += "," + self.__windowedImpulseResp[n].ToString("F9")
				data[3] += "," + self.__windowedStepResp[n].ToString("F9") 
		

			data[1] += "\n\n"
			data[2] += "\n\n"

			System.IO.File.WriteAllLines(saveFileDialog.FileName, data)
			MessageBox.Show("Data written to file!", "Export Time Domain Data", MessageBoxButtons.OK, MessageBoxIcon.Information)


	def __exportFreqDomainDataToolStripMenuItem_Click(self):

		saveFileDialog.Filter = "Text File|*.txt"
		saveFileDialog.Title  = "Export Freq Domain Data"
		saveFileDialog.ShowDialog()

		if (saveFileDialog.FileName != ""):

			# string[] data = new string[4]
			data = [""] * 4
			data[0] = "[FREQUENCY DOMAIN DATA (FREQ/RAW/WINDOWED)] Filter Order: " + self.__NUM_TOTAL_SAMPLES + " Sampling Freq (Hz): " + (1.0 / self.__SAMPLE_TIME_S).ToString("F6") + " Cut-Off Freq Lo (Hz): " + self.__CUTOFF_FREQLO_HZ.ToString("F6") + " Cut-Off Freq Hi (Hz): " + self.__CUTOFF_FREQHI_HZ.ToString("F6") + "\n\n"

			data[1] = self.__freqVecHz[0].ToString("F6")
			data[2] = self.__impRespMag[0].ToString("F9")
			data[3] = self.__winRespMag[0].ToString("F9")
			# for (int n = 1; n < NUM_FREQ_SAMPLES; n++):
			for n in range(self.__NUM_FREQ_SAMPLES):
				data[1] += "," + self.__freqVecHz[n].ToString("F6")
				data[2] += "," + self.__impRespMag[n].ToString("F9")
				data[3] += "," + self.__winRespMag[n].ToString("F9")

			data[1] += "\n\n"
			data[2] += "\n\n"

			System.IO.File.WriteAllLines(saveFileDialog.FileName, data)
			# MessageBox.Show("Data written to file!", "Export Freq Domain Data", MessageBoxButtons.OK, MessageBoxIcon.Information)

if __name__ == '__main__':
	root = tk.Tk()
	myapp = App(root)
	myapp.mainloop()
