#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.pylab import mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk #NavigationToolbar2TkAgg
from math import *
from enum import Enum, unique

# from tkinter import *
import tkinter as tk
import tkinter.messagebox as tkMessageBox
import tkinter.ttk as ttk

from logit import pv


np.seterr(divide='ignore', invalid='ignore')


@unique
class WindowType(Enum):
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


class App(tk.Frame):

	def __init__(self, master):
		super().__init__(master)

		# Initial settings
		
		# Parameters
		self.__SAMPLE_TIME_S = 0.01
		self.__CUTOFF_FREQUENCY_HZ = 20.0
		self.__CUTOFF_FREQUENCY2_HZ = 20.0
		self.__NUM_TOTAL_SAMPLES = 64
		self.__NUM_SHIFT_SAMPLES = 32
		self.__WIN_TYPE = WindowType.Hamming
		self.__FILT_TYPE = FilterType.LowPass
		self.__NUM_FREQ_SAMPLES = 256

		# Time Domain
		# float[] timeVector
		# float[] impulseResponse
		# float[] stepResponse
		# float[] window
		# float[] windowedImpulseResponse
		# float[] windowedStepResponse

		# Frq Domain
		# float[] frequencyVectorHz
		# float[] impRespMag
		# float[] winRespMag
		# float[] winMag

		# cmbWindow.SelectedIndex = 5
		# radViewImpulse.Checked = True
		# cmbDisplay.SelectedIndex = 0

		self.__CreateWindow(master)
		self.__FilterTypeChange()

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
		window.geometry('%dx%d+%d+%d' %(width, hight, cen_x, cen_y))

	def __CreateWindow(self, root):

		# root = tk.Tk()  # 创建tkinter的主窗口
		# root.title("在tkinter中使用matplotlib")

		root.title("Design FIR")  # 设置窗口标题
		# root.geometry("600x600")  # 设置窗口大小 注意：是x 不是*
		self.__CenterWindow(root, 1200, 600)

		# Filter Response in Time Domain
		figureFilterTimeDomain = self.__DrawSin()
		canvasFilterTimeDomain = FigureCanvasTkAgg(figureFilterTimeDomain, root)
		canvasFilterTimeDomain.draw()  #以前的版本使用show()方法，matplotlib 2.2之后不再推荐show（）用draw代替，但是用show不会报错，会显示警告
		# canvasFilterTimeDomain.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		canvasFilterTimeDomain.get_tk_widget().grid(row=0,column=0, columnspan = 4)

		# Window Function in Time Domain
		figureWindowTimeDomain = self.__DrawSin()
		canvasWindowTimeDomain = FigureCanvasTkAgg(figureFilterTimeDomain, root)
		canvasWindowTimeDomain.draw()  #以前的版本使用show()方法，matplotlib 2.2之后不再推荐show（）用draw代替，但是用show不会报错，会显示警告
		# canvasWindowTimeDomain.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		canvasWindowTimeDomain.get_tk_widget().grid(row=0, column=4, columnspan = 6)

		# Filter Frq Response
		figureFilterFrqDomain = self.__DrawSin()
		canvasFilterFrqDomain = FigureCanvasTkAgg(figureFilterFrqDomain, root)
		canvasFilterFrqDomain.draw()  #以前的版本使用show()方法，matplotlib 2.2之后不再推荐show（）用draw代替，但是用show不会报错，会显示警告
		# canvasFilterFrqDomain.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		canvasFilterFrqDomain.get_tk_widget().grid(row=1,column=0, columnspan = 4)

		# Window Function Frq Response
		figureWindowFrqDomain = self.__DrawSin()
		canvasWindowFrqDomain = FigureCanvasTkAgg(figureWindowFrqDomain, root)
		canvasWindowFrqDomain.draw()  #以前的版本使用show()方法，matplotlib 2.2之后不再推荐show（）用draw代替，但是用show不会报错，会显示警告
		# canvasWindowFrqDomain.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		canvasWindowFrqDomain.get_tk_widget().grid(row=1, column=4, columnspan = 6)

		self.__txtSamplingFrq = tk.DoubleVar()
		self.__txtSamplingFrq.set(1000)
		lblSampFrq = tk.Label(text = "Samp Frq(Hz): ")
		lblSampFrq.grid(row=2,column=0)
		enySampFrq = ttk.Entry(textvariable = self.__txtSamplingFrq, width = 10)
		enySampFrq.grid(row=2,column=1)

		self.__txtFilterLength = tk.IntVar()
		self.__txtFilterLength.set(1)
		lblFilterLength = tk.Label(text = "Length of Filter: ")  
		lblFilterLength.grid(row=3,column=0)
		enyFilterLength = ttk.Entry(textvariable = self.__txtFilterLength, width = 10)
		enyFilterLength.grid(row=3,column=1)

		self.__txtShiftSamples = tk.IntVar()
		self.__txtShiftSamples.set(1)
		lblShiftSamples = tk.Label(text = "Shift Samples: ")  
		lblShiftSamples.grid(row=4,column=0)
		enyShiftSamples = ttk.Entry(textvariable = self.__txtShiftSamples, width = 10)
		enyShiftSamples.grid(row=4,column=1)

		self.__txtCutOffFrq = tk.DoubleVar()
		self.__txtCutOffFrq.set(20.0)
		lblCutOffFrq = tk.Label(text = "Low Frq(Hz): ")  
		lblCutOffFrq.grid(row=2,column=2)
		self.__enyCutOffFrq = ttk.Entry(textvariable = self.__txtCutOffFrq, width = 10)
		self.__enyCutOffFrq.grid(row=2,column=3)

		self.__txtCutOffFrqHi = tk.DoubleVar()
		self.__txtCutOffFrqHi.set(20.0)
		lblCutOffFrqHi = tk.Label(text = "High Frq(Hz): ")  
		lblCutOffFrqHi.grid(row=3,column=2)
		self.__enyCutOffFrqHi = ttk.Entry(textvariable = self.__txtCutOffFrqHi, width = 10)
		self.__enyCutOffFrqHi.grid(row=3,column=3)

		self.__cmbWinTyp = tk.StringVar()
		self.__cmbWinTyp.set("0")
		lblcmbWinTyp = tk.Label(text = "Type of Window: ")  
		lblcmbWinTyp.grid(row=4,column=2)
		options = ["Rectangular", "Triangular", "Welch", "Sine", "Hann", "Hamming", "Blackman", "Nuttall", "BlackmanNuttall", "BlackmanHarris", "FlatTop"]
		cmbcmbWinTyp = ttk.Combobox(textvariable = self.__cmbWinTyp, values = options, state='readonly', width = 10)
		cmbcmbWinTyp.current(0)
		cmbcmbWinTyp.bind('<<ComboboxSelected>>', self.__cmbWindow_SelectedIndexChanged)
		cmbcmbWinTyp.grid(row=4,column=3)

		# radBS radBP radHP radLP
		self.__radFiltTyp = tk.IntVar()
		self.__radFiltTyp.set(0)	# 初始化变量值
		lblFiltTyp = tk.Label(text = "Filter Type: ")  
		lblFiltTyp.grid(row=2,column=4, columnspan = 2)
		radio = ttk.Radiobutton(variable = self.__radFiltTyp, value = 0, text ='Low Pass', command = self.__radPS_CheckedChanged)
		radio.grid(row = 3,column = 4, sticky = tk.W)
		radio = ttk.Radiobutton(variable = self.__radFiltTyp, value = 1, text ='High Pass', command = self.__radPS_CheckedChanged)
		radio.grid(row=3,column=5, sticky = tk.W)
		radio = ttk.Radiobutton(variable = self.__radFiltTyp, value = 2, text ='Band Pass', command = self.__radPS_CheckedChanged)
		radio.grid(row=4,column=4, sticky = tk.W)
		radio = ttk.Radiobutton(variable = self.__radFiltTyp, value = 3, text ='Band Stop', command = self.__radPS_CheckedChanged)
		radio.grid(row=4,column=5, sticky = tk.W)

		# radViewStep radViewImpulse
		self.__radRespPltTyp = tk.IntVar()
		self.__radRespPltTyp.set(0)	# 初始化变量值
		lblRespPltTyp = tk.Label(text = "Response Plot: ")  
		lblRespPltTyp.grid(row=2,column=6)
		radio = ttk.Radiobutton(variable = self.__radRespPltTyp, value = 0, text ='Impulse', command = self.__radViewImpulse_CheckedChanged)
		radio.grid(row=2,column=7, sticky = tk.W)
		radio = ttk.Radiobutton(variable = self.__radRespPltTyp, value = 1, text ='Step', command = self.__radViewImpulse_CheckedChanged)
		radio.grid(row=3,column=7, sticky = tk.W)

		# cmbDisplay
		self.__cmbDisplay = tk.StringVar()
		self.__cmbDisplay.set("0")
		lblcmbDisplay = tk.Label(text = "Display: ")  
		lblcmbDisplay.grid(row=4,column=6)
		options = ["Both", "Raw Only", "Win Only"]
		cmbcmbDisplay = ttk.Combobox(textvariable = self.__cmbDisplay, values = options, state='readonly', width = 10)
		cmbcmbDisplay.current(0)
		cmbcmbDisplay.bind('<<ComboboxSelected>>', self.__cmbDisplay_SelectedIndexChanged)
		cmbcmbDisplay.grid(row=4,column=7)		

		# fileToolStripMenuItem
		
		# exportCoefficientsToolStripMenuItem
		
		# exportTimeDomainDataToolStripMenuItem
		
		# exportFrqDomainDataToolStripMenuItem
		
		# infoToolStripMenuItem

		btnDesignFilter = ttk.Button(master = root, text = "Design Filter", command = self.__btnDesignFilter_Click)
		btnDesignFilter.grid(row = 2,column = 8, rowspan = 3, columnspan = 2, sticky = tk.N + tk.S + tk.W + tk.E)

	def __DesignFilter(self):
		timeVector = self.__ComputeTimeVector()
		window = self.__ComputeWindow()
		impulseResponse, stepResponse = self.__ComputeResponses()
		windowedImpulseResponse, windowedStepResponse = self.__ComputeWindowedResponses(impulseResponse, window)

		frequencyVectorHz = self.__ComputeFrqVector()
		impRespMag, winRespMag = self.__ComputeRespBode(impulseResponse, windowedImpulseResponse, frequencyVectorHz)
		winMag = self.__ComputeWindowDFT(window, frequencyVectorHz)	

		# self.__UpdateCharts(timeVector, impulseResponse, windowedImpulseResponse, stepResponse, windowedStepResponse, window, frequencyVectorHz, winRespMag)

	def __UpdateCharts(self, timeVector, impulseResponse, windowedImpulseResponse, stepResponse, windowedStepResponse, window, frequencyVectorHz, winRespMag):

		chrtFilterTimeDomain.Series[0].Points.DataBindXY(timeVector, impulseResponse)
		chrtFilterTimeDomain.Series[1].Points.DataBindXY(timeVector, windowedImpulseResponse)
		chrtFilterTimeDomain.Series[2].Points.DataBindXY(timeVector, stepResponse)
		chrtFilterTimeDomain.Series[3].Points.DataBindXY(timeVector, windowedStepResponse)
		chrtFilterTimeDomain.Series[2].Enabled = False
		chrtFilterTimeDomain.Series[3].Enabled = False

		chrtFilterTimeDomain.ChartAreas[0].RecalculateAxesScale()
		chrtFilterTimeDomain.Update()

		chrtWindowTimeDomain.Series[0].Points.DataBindXY(timeVector, window)
		chrtWindowTimeDomain.ChartAreas[0].RecalculateAxesScale()

		chrtFilterFrqDomain.Series[0].Points.DataBindXY(frequencyVectorHz, impRespMag)
		chrtFilterFrqDomain.Series[1].Points.DataBindXY(frequencyVectorHz, winRespMag)
		chrtFilterFrqDomain.ChartAreas[0].AxisX.Interval = 10.0 * (0.5 / self.__SAMPLE_TIME_S) / (self.__NUM_FREQ_SAMPLES - 1.0)
		chrtFilterFrqDomain.ChartAreas[0].RecalculateAxesScale()

		chrtWindowFrqDomain.Series[0].Points.DataBindXY(frequencyVectorHz, winMag)
		chrtWindowFrqDomain.ChartAreas[0].AxisX.Interval = 10.0 * (0.5 / self.__SAMPLE_TIME_S) / (self.__NUM_FREQ_SAMPLES - 1.0)
		chrtWindowFrqDomain.ChartAreas[0].RecalculateAxesScale()


	# Time Domain Functions
	def __ComputeTimeVector(self):
		# timeVector = new float[self.__NUM_TOTAL_SAMPLES]
		timeVector = [0.0] * self.__NUM_TOTAL_SAMPLES
		# for (int n = 0 n < self.__NUM_TOTAL_SAMPLES n++):
		for n in range(self.__NUM_TOTAL_SAMPLES):
			timeVector[n] = n * self.__SAMPLE_TIME_S
	
		return timeVector
	
	def __ComputeResponses(self):

		# impulseResponse = new float[self.__NUM_TOTAL_SAMPLES]
		impulseResponse = [0.0] * self.__NUM_TOTAL_SAMPLES
		# stepResponse = new float[self.__NUM_TOTAL_SAMPLES]
		stepResponse = [0.0] * self.__NUM_TOTAL_SAMPLES

		# for (int n = 0 n < self.__NUM_TOTAL_SAMPLES n++)
		for n in range(self.__NUM_TOTAL_SAMPLES):
			if (n != self.__NUM_SHIFT_SAMPLES):
				if self.__FILT_TYPE == FilterType.LowPass:
					impulseResponse[n] = math.sin(2.0 * math.pi * self.__CUTOFF_FREQUENCY_HZ * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)) / (math.pi * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES))
				elif self.__FILT_TYPE == FilterType.HighPass:							
					impulseResponse[n] = (math.sin(math.pi * (n - self.__NUM_SHIFT_SAMPLES)) - math.sin(2.0 * math.pi * self.__CUTOFF_FREQUENCY_HZ * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES))) / (math.pi * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES))
				elif self.__FILT_TYPE == FilterType.BandPass:
					impulseResponse[n] = (math.sin(2.0 * math.pi * self.__CUTOFF_FREQUENCY2_HZ * self.__SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES)) - math.sin(2.0 * math.pi * self.__CUTOFF_FREQUENCY_HZ * self.__SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES))) / (math.pi * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES))
				elif self.__FILT_TYPE == FilterType.BandStop:
					impulseResponse[n] = (math.sin(2.0 * math.pi * self.__CUTOFF_FREQUENCY_HZ * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)) - math.sin(2.0 * math.pi * self.__CUTOFF_FREQUENCY2_HZ * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES)) + math.sin(math.pi * (n - self.__NUM_SHIFT_SAMPLES))) / (math.pi * self.__SAMPLE_TIME_S * (n - self.__NUM_SHIFT_SAMPLES))
			
			else: # Avoid divide-by-zero, limit is 2*fc
				if self.__FILT_TYPE == FilterType.LowPass:
					impulseResponse[n] = 2.0 * self.__CUTOFF_FREQUENCY_HZ
					break
				elif self.__FILT_TYPE == FilterType.HighPass:
					impulseResponse[n] = 1.0 / self.__SAMPLE_TIME_S - 2.0 * self.__CUTOFF_FREQUENCY_HZ
					break
				elif self.__FILT_TYPE == FilterType.BandPass:
					impulseResponse[n] = 2.0 * self.__CUTOFF_FREQUENCY2_HZ - 2.0 * self.__CUTOFF_FREQUENCY_HZ
					break
				elif self.__FILT_TYPE == FilterType.BandStop:
					impulseResponse[n] = 2.0 * self.__CUTOFF_FREQUENCY_HZ - 2.0 * self.__CUTOFF_FREQUENCY2_HZ + 1.0 / self.__SAMPLE_TIME_S
					break

		# Normalise by DC gain to achieve 0dB gain at DC and then compute step response
		# for (int n = 0 n < self.__NUM_TOTAL_SAMPLES n++)
		for n in range(self.__NUM_TOTAL_SAMPLES):
			impulseResponse[n] *= self.__SAMPLE_TIME_S

			if (n == 0):
				stepResponse[n] = 0.5 * impulseResponse[n]
			else:		
				stepResponse[n] = stepResponse[n - 1] + 0.5 * (impulseResponse[n] + impulseResponse[n - 1])
		return impulseResponse, stepResponse


	def __ComputeWindow(self):

		# window = new float[self.__NUM_TOTAL_SAMPLES]
		window = [0.0] * self.__NUM_TOTAL_SAMPLES

		# for (int n = 0 n < self.__NUM_TOTAL_SAMPLES n++)
		for n in range(self.__NUM_TOTAL_SAMPLES):
			if self.__WIN_TYPE == "Rectangular":
				window[n] = 1.0
			elif self.__WIN_TYPE == "Triangular":
				window[n] = 1.0 - Abs((n - 0.5 * self.__NUM_TOTAL_SAMPLES) / (0.5 * self.__NUM_TOTAL_SAMPLES))
			elif self.__WIN_TYPE == "Welch":
				window[n] = 1.0 - Pow((n - 0.5 * self.__NUM_TOTAL_SAMPLES) / (0.5 * self.__NUM_TOTAL_SAMPLES), 2.0)
			elif self.__WIN_TYPE == "Sine":
				window[n] = math.sin(math.pi * n / (self.__NUM_TOTAL_SAMPLES))
			elif self.__WIN_TYPE == "Hann":
				window[n] = 0.5 * (1 - math.math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)))
			elif self.__WIN_TYPE == "Hamming":
				window[n] = (25.0 / 46.0) - (21.0 / 46.0) * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES))
			elif self.__WIN_TYPE == "Blackman":
				window[n] = 0.42 - 0.5 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.08 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES))
			elif self.__WIN_TYPE == "Nuttall":
				window[n] = 0.355768 - 0.487396 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.144232 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) - 0.012604 * math.cos(6.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES))
			elif self.__WIN_TYPE == "BlackmanNuttall":
				window[n] = 0.3635819 - 0.4891775 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.1365995 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) - 0.0106411 * math.cos(6.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES))
			elif self.__WIN_TYPE == "BlackmanHarris":
				window[n] = 0.35875 - 0.48829 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.14128 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) - 0.01168 * math.cos(6.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES))
			elif self.__WIN_TYPE == "FlatTop":
				window[n] = 0.21557895 - 0.41663158 * math.cos(2.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.277263158 * math.cos(4.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) - 0.083578947 * math.cos(6.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES)) + 0.006947368 * math.cos(8.0 * math.pi * n / (self.__NUM_TOTAL_SAMPLES))
			else:
				window[n] = 1.0
		return window

	def __ComputeWindowedResponses(self, impulseResponse, window):

		# windowedImpulseResponse = new float[self.__NUM_TOTAL_SAMPLES]
		windowedImpulseResponse = [0.0] * self.__NUM_TOTAL_SAMPLES
		# windowedStepResponse = new float[self.__NUM_TOTAL_SAMPLES]
		windowedStepResponse = [0.0] * self.__NUM_TOTAL_SAMPLES

		# for (int n = 0; n < self.__NUM_TOTAL_SAMPLES; n++):
		for n in range(self.__NUM_TOTAL_SAMPLES):
			windowedImpulseResponse[n] = impulseResponse[n] * window[n]

			if (n == 0):
				windowedStepResponse[n] = 0.5 * windowedStepResponse[n]
			else:		
				windowedStepResponse[n] = windowedStepResponse[n - 1] + 0.5 * (windowedImpulseResponse[n] + windowedImpulseResponse[n - 1])

		return windowedImpulseResponse, windowedStepResponse


	# Frq Domain Functions
	def __ComputeFrqVector(self):

		# frequencyVectorHz = new float[self.__NUM_FREQ_SAMPLES]
		frequencyVectorHz = [0.0] * self.__NUM_FREQ_SAMPLES

		df = (0.5 / self.__SAMPLE_TIME_S) / (self.__NUM_FREQ_SAMPLES - 1.0)

		# for (int n = 0; n < self.__NUM_FREQ_SAMPLES; n++):
		for n in range(self.__NUM_TOTAL_SAMPLES):
			frequencyVectorHz[n] = n * df

		return frequencyVectorHz

	def __ComputeRespBode(self, impulseResponse, windowedImpulseResponse, frequencyVectorHz):

		# impRespMag = new float[self.__NUM_FREQ_SAMPLES]
		impRespMag = [0.0] * self.__NUM_FREQ_SAMPLES
		# winRespMag = new float[self.__NUM_FREQ_SAMPLES]
		winRespMag = [0.0] * self.__NUM_FREQ_SAMPLES

		# for (int fIndex = 0; fIndex < NUM_FREQ_SAMPLES; fIndex++):
		for fIndex in range(self.__NUM_FREQ_SAMPLES):
			re = 0.0
			im = 0.0
			reWin = 0.0
			imWin = 0.0

			# for (int n = 0; n < self.__NUM_TOTAL_SAMPLES; n++)
			for n in range(self.__NUM_TOTAL_SAMPLES):
				re = re + impulseResponse[n] * math.cos(2.0 * math.pi * frequencyVectorHz[fIndex] * n * self.__SAMPLE_TIME_S)
				im = im - impulseResponse[n] * math.sin(2.0 * math.pi * frequencyVectorHz[fIndex] * n * self.__SAMPLE_TIME_S)
				reWin = reWin + windowedImpulseResponse[n] * math.cos(2.0 * math.pi * frequencyVectorHz[fIndex] * n * self.__SAMPLE_TIME_S)
				imWin = imWin - windowedImpulseResponse[n] * math.sin(2.0 * math.pi * frequencyVectorHz[fIndex] * n * self.__SAMPLE_TIME_S)

			impRespMag[fIndex] = 10.0 * math.log10(re * re + im * im)
			pv(reWin * reWin + imWin * imWin)
			winRespMag[fIndex] = 10.0 * math.log10(reWin * reWin + imWin * imWin)

		return impRespMag, winRespMag

	def __GetGainAtCutOff(self) -> float: 
		re = 0.0
		im = 0.0

		# for (int n = 0; n < self.__NUM_TOTAL_SAMPLES; n++):
		for n in range(self.__NUM_TOTAL_SAMPLES):
			re = re + impulseResponse[n] * math.cos(2.0 * math.pi * self.__CUTOFF_FREQUENCY_HZ * n * self.__SAMPLE_TIME_S)
			im = im - impulseResponse[n] * math.sin(2.0 * math.pi * self.__CUTOFF_FREQUENCY_HZ * n * self.__SAMPLE_TIME_S)

		return (10.0 * Log10(re * re + im * im))


	def __ComputeWindowDFT(self, window, frequencyVectorHz):

		# winMag = new float[self.__NUM_FREQ_SAMPLES]
		winMag = [0.0] * self.__NUM_FREQ_SAMPLES

		# for (int fIndex = 0; fIndex < NUM_FREQ_SAMPLES; fIndex++):
		for fIndex in range(self.__NUM_FREQ_SAMPLES):

			re = 0.0
			im = 0.0

			# for (int n = 0; n < NUM_TOTAL_SAMPLES; n++):
			for n in range(self.__NUM_TOTAL_SAMPLES):
				re = re + window[n] * math.cos(2.0 * math.pi * frequencyVectorHz[fIndex] * n * self.__SAMPLE_TIME_S)
				im = im - window[n] * math.sin(2.0 * math.pi * frequencyVectorHz[fIndex] * n * self.__SAMPLE_TIME_S)

			winMag[fIndex] = 10.0 * math.log10(re * re + im * im)

		return winMag

	def __UpdatePlotSettings(self):

		if (radViewImpulse.Checked):	
			if (cmbDisplay.SelectedIndex == 0):		
				chrtFilterTimeDomain.Series[0].configure(state='normal')
				chrtFilterTimeDomain.Series[1].configure(state='normal')
			elif (cmbDisplay.SelectedIndex == 1):		
				chrtFilterTimeDomain.Series[0].configure(state='normal')
				chrtFilterTimeDomain.Series[1].Enabled = False		
			elif (cmbDisplay.SelectedIndex == 2):		
				chrtFilterTimeDomain.Series[0].Enabled = False
				chrtFilterTimeDomain.Series[1].configure(state='normal')
			chrtFilterTimeDomain.Series[2].Enabled = False
			chrtFilterTimeDomain.Series[3].Enabled = False

		else:
			if (cmbDisplay.SelectedIndex == 0):		
				chrtFilterTimeDomain.Series[2].configure(state='normal')
				chrtFilterTimeDomain.Series[3].configure(state='normal')		
			elif (cmbDisplay.SelectedIndex == 1):		
				chrtFilterTimeDomain.Series[2].configure(state='normal')
				chrtFilterTimeDomain.Series[3].Enabled = False		
			elif (cmbDisplay.SelectedIndex == 2):		
				chrtFilterTimeDomain.Series[2].Enabled = False
				chrtFilterTimeDomain.Series[3].configure(state='normal')
			chrtFilterTimeDomain.Series[0].Enabled = False
			chrtFilterTimeDomain.Series[1].Enabled = False

		if (cmbDisplay.SelectedIndex == 0):	
			chrtFilterFrqDomain.Series[0].configure(state='normal')
			chrtFilterFrqDomain.Series[1].configure(state='normal')	
		elif (cmbDisplay.SelectedIndex == 1):	
			chrtFilterFrqDomain.Series[0].configure(state='normal')
			chrtFilterFrqDomain.Series[1].Enabled = False	
		elif (cmbDisplay.SelectedIndex == 2):	
			chrtFilterFrqDomain.Series[0].Enabled = False
			chrtFilterFrqDomain.Series[1].configure(state='normal')

		chrtFilterTimeDomain.ChartAreas[0].RecalculateAxesScale()

	def create_matplotlib(self):
		#创建绘图对象f
		f=plt.figure(num=2,figsize=(16,12),dpi=80,facecolor="pink",edgecolor='green',frameon=True)
		#创建一副子图
		fig1=plt.subplot(1,1,1)
 
		x=np.arange(0,2*np.pi,0.1)
		y1=np.sin(x)
		y2=np.cos(x)
 
		line1,=fig1.plot(x,y1,color='red',linewidth=3,linestyle='--')	#画第一条线
		line2,=fig1.plot(x,y2) 
		plt.setp(line2,color='black',linewidth=8,linestyle='-',alpha=0.3)#华第二条线
 
		fig1.set_title("这是第一幅图",loc='center',pad=20,fontsize='xx-large',color='red')	#设置标题
		line1.set_label("正弦曲线")														   #确定图例
		fig1.legend(['正弦','余弦'],loc='upper left',facecolor='green',frameon=True,shadow=True,framealpha=0.5,fontsize='xx-large')
 
		fig1.set_xlabel('横坐标')															 #确定坐标轴标题
		fig1.set_ylabel("纵坐标")
		fig1.set_yticks([-1,-1/2,0,1/2,1])												   #设置坐标轴刻度
		fig1.grid(which='major',axis='x',color='r', linestyle='-', linewidth=2)			  #设置网格
		
		return f

	def __DrawSin(self):
		f = plt.Figure(figsize=(6.0, 2.5), dpi=100)
		a = f.add_subplot(111)  # 添加子图:1行1列第1个

		# 生成用于绘sin图的数据
		x = np.arange(0, 3, 0.01)
		y = np.sin(2 * np.pi * x)

		# 在前面得到的子图上绘图
		a.plot(x, y)

		return f;

	def __radViewImpulse_CheckedChanged(self):
		self.__UpdatePlotSettings()

	def __cmbDisplay_SelectedIndexChanged(self):
		self.__UpdatePlotSettings()

	def __btnDesignFilter_Click(self):

		# self.__SAMPLE_TIME_S = 1.0 / float(txtSamplingFrq.get())
		sampFrq = float(self.__txtSamplingFrq.get())
		if sampFrq == 0:
			self.__ShowErr("Sampling frequency cannot be zero")
			return
		self.__SAMPLE_TIME_S = 1.0 / sampFrq
		self.__NUM_TOTAL_SAMPLES = self.__txtFilterLength.get()
		self.__NUM_SHIFT_SAMPLES = self.__txtShiftSamples.get()
		self.__FILT_TYPE = FilterType(self.__radFiltTyp.get())
		pv(self.__FILT_TYPE)
		self.__WIN_TYPE = self.__cmbWinTyp.get()
		pv(self.__WIN_TYPE)

		if (self.__FILT_TYPE == FilterType.LowPass):
			self.__CUTOFF_FREQUENCY_HZ = float(self.__txtCutOffFrqHigh.get())
		elif (self.__FILT_TYPE == FilterType.HighPass):
			self.__CUTOFF_FREQUENCY_HZ = float(self.__txtCutOffFrq.get())
		elif (self.__FILT_TYPE == FilterType.BandPass):
			self.__CUTOFF_FREQUENCY_HZ = float(self.__txtCutOffFrq.get())
			self.__CUTOFF_FREQUENCY2_HZ = float(self.__txtCutOffFrqHigh.get())
		elif (self.__FILT_TYPE == FilterType.BandStop):
			self.__CUTOFF_FREQUENCY_HZ = float(self.__txtCutOffFrq.get())
			self.__CUTOFF_FREQUENCY2_HZ = float(self.__txtCutOffFrqHigh.get())

		if (self.__SAMPLE_TIME_S < 0.0):
			self.__ShowErr("Sampling frequency cannot be negative.")
			return

		if (self.__CUTOFF_FREQUENCY_HZ >= 0.5 / self.__SAMPLE_TIME_S or self.__CUTOFF_FREQUENCY2_HZ >= 0.5 / self.__SAMPLE_TIME_S):	
			self.__ShowErr("Cut-off frequency has to be less than the Nyquist frequency (i.e. sampling frequency / 2).")
			return

		if (self.__NUM_TOTAL_SAMPLES < 0 or self.__NUM_SHIFT_SAMPLES < 0):	
			self.__ShowErr("Total number of samples and sample shift number both need to be integers, greater than zero.")
			return

		self.__DesignFilter()

		# self.__UpdatePlotSettings()

	def __cmbWindow_SelectedIndexChanged(self, p):

		if (cmbWindow.SelectedIndex >= 0):
	
			# self.__WIN_TYPE = (WindowType)Convert.ToInt32(cmbWindow.SelectedIndex)

			self.__ComputeFrqVector()
			self.__ComputeWindow()
			self.__ComputeWindowDFT()

			chrtWindowTimeDomain.Series[0].Points.DataBindXY(timeVector, window)

			chrtWindowFrqDomain.Series[0].Points.DataBindXY(frequencyVectorHz, winMag)

			chrtWindowTimeDomain.ChartAreas[0].RecalculateAxesScale()
			chrtWindowFrqDomain.ChartAreas[0].RecalculateAxesScale()
	
	def __infoToolStripMenuItem_Click(self):

		MessageBox.Show("FIR Filter Designer\nWritten by Philip M. Salmony\n29 November 2019\nphilsal.co.uk", "Info", MessageBoxButtons.OK, MessageBoxIcon.Information)


	def __FilterTypeChange(self):
		stFiltTyp = FilterType(self.__radFiltTyp.get())
		pv(stFiltTyp)

		if (stFiltTyp == FilterType.LowPass):	
			self.__enyCutOffFrq.configure(state='disabled')
			self.__enyCutOffFrqHi.configure(state='normal')
		elif (stFiltTyp == FilterType.HighPass):	
			self.__enyCutOffFrq.configure(state='normal')
			self.__enyCutOffFrqHi.configure(state='disabled')
		elif (stFiltTyp == FilterType.BandPass):	
			self.__enyCutOffFrq.configure(state='normal')
			self.__enyCutOffFrqHi.configure(state='normal')
		elif (stFiltTyp == FilterType.BandStop):
			self.__enyCutOffFrq.configure(state='normal')
			self.__enyCutOffFrqHi.configure(state='normal')

	def __radPS_CheckedChanged(self):
		self.__FilterTypeChange()

	def __exportCoefficientsToolStripMenuItem_Click(self):

		saveFileDialog.Filter = "Text File|*.txt"
		saveFileDialog.Title  = "Export Filter Coefficients"
		saveFileDialog.ShowDialog()

		if (saveFileDialog.FileName != ""):
	
			# string[] data = new string[3]
			data = [""] * 3
			data[0] = "Filter Order: " + self.__NUM_TOTAL_SAMPLES + " Sampling Frq (Hz): " + (1.0 / SAMPLE_TIME_S).ToString("F6") + " Cut-Off Frq Lo (Hz): " + self.__CUTOFF_FREQUENCY_HZ.ToString("F6") + " Cut-Off Frq Hi (Hz): " + self.__CUTOFF_FREQUENCY2_HZ.ToString("F6") + "\n\n"

			data[1] = windowedImpulseResponse[0].ToString("F7")
			data[2] = "float coeff[] = {" + windowedImpulseResponse[0].ToString("F7") + "f"
			# for (int n = 1 n < NUM_TOTAL_SAMPLES n++)
			for n in range(self.__NUM_TOTAL_SAMPLES):
				data[1] += "," + windowedImpulseResponse[n].ToString("F9")
				data[2] += "," + windowedImpulseResponse[n].ToString("F7") + "f"
		
			data[1] += "\n\n"
			data[2] += "}"
			
			System.IO.File.WriteAllLines(saveFileDialog.FileName, data)
			MessageBox.Show("Coefficients written to file!", "Export Coefficients", MessageBoxButtons.OK, MessageBoxIcon.Information)
	

	def __exportTimeDomainDataToolStripMenuItem_Click(self):

		saveFileDialog.Filter = "Text File|*.txt"
		saveFileDialog.Title  = "Export Time Domain Data"
		saveFileDialog.ShowDialog()

		if (saveFileDialog.FileName != ""):
	
			# string[] data = new string[4]
			data = [""] * 4
			data[0] = "[TIME DOMAIN DATA (TIME/IMPULSE/STEP)] Filter Order: " + self.__NUM_TOTAL_SAMPLES + " Sampling Frq (Hz): " + (1.0 / self.__SAMPLE_TIME_S).ToString("F6") + " Cut-Off Frq Lo (Hz): " + self.__CUTOFF_FREQUENCY_HZ.ToString("F6") + " Cut-Off Frq Hi (Hz): " + self.__CUTOFF_FREQUENCY2_HZ.ToString("F6") + "\n\n"

			data[1] = timeVector[0].ToString("F6")
			data[2] = windowedImpulseResponse[0].ToString("F9")
			data[3] = windowedStepResponse[0].ToString("F9")

			# for (int n = 1; n < NUM_TOTAL_SAMPLES; n++):
			for n in range(self.__NUM_TOTAL_SAMPLES):
				data[1] += "," + timeVector[n].ToString("F6")
				data[2] += "," + windowedImpulseResponse[n].ToString("F9")
				data[3] += "," + windowedStepResponse[n].ToString("F9") 
		

			data[1] += "\n\n"
			data[2] += "\n\n"

			System.IO.File.WriteAllLines(saveFileDialog.FileName, data)
			MessageBox.Show("Data written to file!", "Export Time Domain Data", MessageBoxButtons.OK, MessageBoxIcon.Information)


	def __exportFrqDomainDataToolStripMenuItem_Click(self):

		saveFileDialog.Filter = "Text File|*.txt"
		saveFileDialog.Title  = "Export Frq Domain Data"
		saveFileDialog.ShowDialog()

		if (saveFileDialog.FileName != ""):

			# string[] data = new string[4]
			data = [""] * 4
			data[0] = "[FREQUENCY DOMAIN DATA (FREQ/RAW/WINDOWED)] Filter Order: " + self.__NUM_TOTAL_SAMPLES + " Sampling Frq (Hz): " + (1.0 / self.__SAMPLE_TIME_S).ToString("F6") + " Cut-Off Frq Lo (Hz): " + self.__CUTOFF_FREQUENCY_HZ.ToString("F6") + " Cut-Off Frq Hi (Hz): " + self.__CUTOFF_FREQUENCY2_HZ.ToString("F6") + "\n\n"

			data[1] = frequencyVectorHz[0].ToString("F6")
			data[2] = impRespMag[0].ToString("F9")
			data[3] = winRespMag[0].ToString("F9")
			# for (int n = 1; n < NUM_FREQ_SAMPLES; n++):
			for n in range(self.__NUM_FREQ_SAMPLES):
				data[1] += "," + frequencyVectorHz[n].ToString("F6")
				data[2] += "," + impRespMag[n].ToString("F9")
				data[3] += "," + winRespMag[n].ToString("F9")

			data[1] += "\n\n"
			data[2] += "\n\n"

			System.IO.File.WriteAllLines(saveFileDialog.FileName, data)
			# MessageBox.Show("Data written to file!", "Export Frq Domain Data", MessageBoxButtons.OK, MessageBoxIcon.Information)

if __name__ == '__main__':
	root = tk.Tk()
	myapp = App(root)
	myapp.mainloop()