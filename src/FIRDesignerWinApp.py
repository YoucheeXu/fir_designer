#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import sys
import os
import math
from enum import IntEnum, unique, Enum

import numpy as np

from pyutilities.logit import pv
from pyutilities.matplot import LineData
from pyutilities.tkwin import tkWin


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

class FIRDesignerApp(tkWin):

    def __init__(self):
        super().__init__()

        # Initial settings
        self.__numFreqSample = 256

        # cmbWin.SelectedIndex = 5
        # radViewImpulse.Checked = True
        # cmbDisplay.SelectedIndex = 0

    def init_paras(self):
        self._varSampleTime.set("0.01")
        self._varLowFrequency.set("10.0")
        self._varHighFrequency.set("20.0")
        self._varFilterLength.set("64")
        self._varShiftSamples.set("32")
        self._varWindowType.set("Hamming")
        self._varFilterType.set(0)	# FilrType.LowPass
        self._varResponsePlot.set(0)	# Impulse
        self._varDisplayType.set("Both")		

        # Time Domain
        numTotalSample = int(self._varFilterLength.get())
        self.__timeVec = [0.0] * numTotalSample
        self.__impulseResp = [0.0] * numTotalSample
        self.__stepResp = [0.0] * numTotalSample
        self.__window = [0.0] * numTotalSample
        self.__windowedImpResp = [0.0] * numTotalSample
        self.__windowedStepResp = [0.0] * numTotalSample

        # Frequency Domain
        numFreqSample = self.__numFreqSample
        self.__freqVec = [0.0] * numFreqSample
        self.__impRespMag = [0.0] * numFreqSample
        self.__winRespMag = [0.0] * numFreqSample
        self.__winMag = [0.0] * numFreqSample

    def go(self):
        self._varFilterTypeChanged()
        self.__DesignFilter()
        self._frmApp.mainloop()

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

        self._FilterResponseinTimeDomain.set_xData(self.__timeVec)
        yLine = LineData(self.__impulseResp, "Impulse Response")
        self._FilterResponseinTimeDomain.add_line(yLine, 0)
        yLine = LineData(self.__windowedImpResp, "Windowed Impulse Response")
        self._FilterResponseinTimeDomain.add_line(yLine, 1)
        yLine = LineData(self.__stepResp, "Step Response", visible = False)
        self._FilterResponseinTimeDomain.add_line(yLine, 2)
        yLine = LineData(self.__windowedStepResp, "Windowed Step Response", visible = False)
        self._FilterResponseinTimeDomain.add_line(yLine, 3)
        self._FilterResponseinTimeDomain.draw()

        self._WindowFunctioninTimeDomain.set_xData(self.__timeVec)
        yLine = LineData(self.__window)
        self._WindowFunctioninTimeDomain.add_line(yLine, 0)
        self._WindowFunctioninTimeDomain.draw()

        self._FilterFrequencyResponse.set_xData(self.__freqVec)
        yLine = LineData(self.__impRespMag)
        self._FilterFrequencyResponse.add_line(yLine, 0)
        yLine = LineData(self.__winRespMag)
        self._FilterFrequencyResponse.add_line(yLine, 1)
        self._FilterFrequencyResponse.draw()

        self._WindowFunctionFrequencyResponse.set_xData(self.__freqVec)
        yLine = LineData(self.__winMag)
        self._WindowFunctionFrequencyResponse.add_line(yLine, 0)
        self._WindowFunctionFrequencyResponse.draw()

    # Time Domain Functions
    def __ComputeTimeVec(self):
        numTotalSample = int(self._varFilterLength.get())
        # self.__timeVec = new float[self.__numTotalSample]
        self.__timeVec = [0.0] * numTotalSample

        tiSample = float(self._varSampleTime.get())
        # for (int n = 0 n < self.__numTotalSample n++):
        for n in range(numTotalSample):
            self.__timeVec[n] = n * tiSample

        # pv(self.__timeVec)

    def __ComputeResps(self):
        numTotalSample = int(self._varFilterLength.get())
        # self.__impulseResp = new float[self.__numTotalSample]
        self.__impulseResp = [0.0] * numTotalSample
        # self.__stepResp = new float[self.__numTotalSample]
        self.__stepResp = [0.0] * numTotalSample

        numShftSample = int(self._varShiftSamples.get())
        stFiltTyp = FilrType(self._varFilterType.get())
        tiSample = float(self._varSampleTime.get())
        frqCutOffLo = float(self._varLowFrequency.get())
        frqCutOffHi = float(self._varHighFrequency.get())
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
        numTotalSample = int(self._varFilterLength.get())
        # self.__window = new float[self.__numTotalSample]
        self.__window = [0.0] * numTotalSample

        stWinTyp = self._varWindowType.get()

        # for (int n = 0 n < numTotalSample n++)
        for n in range(numTotalSample):
            if stWinTyp == "Rectangular":
                self.__window[n] = 1.0
            elif stWinTyp == "Triangular":
                self.__window[n] = 1.0 - abs((n - 0.5 * numTotalSample) / (0.5 * numTotalSample))
            elif stWinTyp == "Welch":
                self.__window[n] = 1.0 - math.pow((n - 0.5 * numTotalSample) / (0.5 * numTotalSample), 2.0)
            elif stWinTyp == "Sine":
                self.__window[n] = math.sin(math.pi * n / (numTotalSample))
            elif stWinTyp == "Hann":
                self.__window[n] = 0.5 * (1 - math.cos(2.0 * math.pi * n / (numTotalSample)))
            elif stWinTyp == "Hamming":
                self.__window[n] = (25.0 / 46.0) - (21.0 / 46.0) * math.cos(2.0 * math.pi * n / (numTotalSample))
            elif stWinTyp == "Blackman":
                self.__window[n] = 0.42 - 0.5 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.08 * math.cos(4.0 * math.pi * n / (numTotalSample))
            elif stWinTyp == "Nuttall":
                self.__window[n] = 0.355768 - 0.487396 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.144232 * math.cos(4.0 * math.pi * n / (numTotalSample)) - 0.012604 * math.cos(6.0 * math.pi * n / (numTotalSample))
            elif stWinTyp == "BlackmanNuttall":
                self.__window[n] = 0.3635819 - 0.4891775 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.1365995 * math.cos(4.0 * math.pi * n / (numTotalSample)) - 0.0106411 * math.cos(6.0 * math.pi * n / (numTotalSample))
            elif stWinTyp == "BlackmanHarris":
                self.__window[n] = 0.35875 - 0.48829 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.14128 * math.cos(4.0 * math.pi * n / (numTotalSample)) - 0.01168 * math.cos(6.0 * math.pi * n / (numTotalSample))
            elif stWinTyp == "FlatTop":
                self.__window[n] = 0.21557895 - 0.41663158 * math.cos(2.0 * math.pi * n / (numTotalSample)) + 0.277263158 * math.cos(4.0 * math.pi * n / (numTotalSample)) - 0.083578947 * math.cos(6.0 * math.pi * n / (numTotalSample)) + 0.006947368 * math.cos(8.0 * math.pi * n / (numTotalSample))
            else:
                self.__window[n] = 1.0

        # pv(self.__window)

    def __ComputeWindowedResps(self):
        numTotalSample = int(self._varFilterLength.get())
        # self.__windowedImpResp = new float[self.__numTotalSample]
        self.__windowedImpResp = [0.0] * numTotalSample
        # self.__windowedStepResp = new float[self.__numTotalSample]
        self.__windowedStepResp = [0.0] * numTotalSample

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

        tiSample = float(self._varSampleTime.get())
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

        numTotalSample = int(self._varFilterLength.get())
        tiSample = float(self._varSampleTime.get())
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

        numTotalSample = int(self._varFilterLength.get())
        # for (int n = 0; n < self.__numTotalSample; n++):
        for n in range(numTotalSample):
            re = re + self.__impulseResp[n] * math.cos(2.0 * math.pi * frqCutOffLo * n * tiSample)
            im = im - self.__impulseResp[n] * math.sin(2.0 * math.pi * frqCutOffLo * n * tiSample)

        return (10.0 * math.log10(re * re + im * im))

    def __ComputeWindowDFT(self):
        numFreqSample = self.__numFreqSample
        # winMag = new float[self.__numFreqSample]
        self.__winMag = [0.0] * numFreqSample

        numTotalSample = int(self._varFilterLength.get())
        tiSample = float(self._varSampleTime.get())
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
        stRespPltTyp = RespPltType(self._varResponsePlot.get())
        pv(stRespPltTyp)
        stDisplay = self._varDisplayType.get()
        pv(stDisplay)

        if (stRespPltTyp == RespPltType.Impulse):		# Impulse
            if (stDisplay == "Both"):
                self._FilterResponseinTimeDomain.show_line(0, True)
                self._FilterResponseinTimeDomain.show_line(1, True)
            elif (stDisplay == "Raw Only"):
                self._FilterResponseinTimeDomain.show_line(0, True)
                self._FilterResponseinTimeDomain.show_line(1, False)
            elif (stDisplay == "Win Only"):
                self._FilterResponseinTimeDomain.show_line(0, False)
                self._FilterResponseinTimeDomain.show_line(1, True)				
            self._FilterResponseinTimeDomain.show_line(2, False)
            self._FilterResponseinTimeDomain.show_line(3, False)			

        else:		# Step
            if (stDisplay == "Both"):
                self._FilterResponseinTimeDomain.show_line(2, True)
                self._FilterResponseinTimeDomain.show_line(3, True)
            elif (stDisplay == "Raw Only"):
                self._FilterResponseinTimeDomain.show_line(2, True)
                self._FilterResponseinTimeDomain.show_line(3, False)
            elif (stDisplay == "Win Only"):
                self._FilterResponseinTimeDomain.show_line(2, False)
                self._FilterResponseinTimeDomain.show_line(3, True)
            self._FilterResponseinTimeDomain.show_line(0, False)
            self._FilterResponseinTimeDomain.show_line(1, False)

        if (stDisplay == "Both"):
            self._FilterFrequencyResponse.show_line(0, True)
            self._FilterFrequencyResponse.show_line(1, True)
        elif (stDisplay == "Raw Only"):
            self._FilterFrequencyResponse.show_line(0, True)
            self._FilterFrequencyResponse.show_line(1, False)
        elif (stDisplay == "Win Only"):
            self._FilterFrequencyResponse.show_line(0, False)
            self._FilterFrequencyResponse.show_line(1, True)

        self._FilterResponseinTimeDomain.draw()

        self._FilterFrequencyResponse.draw()

    def _varResponsePlotChanged(self):
        self.__UpdatePlotSettings()

    def _varDisplayTypeChanged(self, event):
        self.__UpdatePlotSettings()

    def _btnDesignFilterClick(self):
        tiSample = float(self._varSampleTime.get())
        pv(tiSample)
        numTotalSample = int(self._varFilterLength.get())
        pv(numTotalSample)
        numShftSample = int(self._varShiftSamples.get())
        pv(numShftSample)
        stWinTyp = self._varWindowType.get()
        pv(stWinTyp)
        stFiltTyp = FilrType(self._varFilterType.get())
        pv(stFiltTyp)

        frqCutOffLo = 0
        frqCutOffHi = 0
        if (stFiltTyp == FilrType.LowPass):
            frqCutOffLo = float(self._varLowFrequency.get())
            pv(frqCutOffLo)
        elif (stFiltTyp == FilrType.HighPass):
            frqCutOffHi = float(self._varHighFrequency.get())
            pv(frqCutOffHi)
        elif (stFiltTyp == FilrType.BandPass):
            frqCutOffLo = float(self._varLowFrequency.get())
            frqCutOffHi = float(self._varHighFrequency.get())
            pv(frqCutOffLo)
            pv(frqCutOffHi)
        elif (stFiltTyp == FilrType.BandStop):
            frqCutOffLo = float(self._varLowFrequency.get())
            frqCutOffHi = float(self._varHighFrequency.get())
            pv(frqCutOffLo)
            pv(frqCutOffHi)

        if (tiSample < 0.0):
            self.ShowErr("Sampling frequency cannot be negative.")
            return

        if (frqCutOffLo >= 0.5 / tiSample or frqCutOffHi >= 0.5 / tiSample):
            self.ShowErr("Cut-off frequency has to be less than the Nyquist frequency (i.e. sampling freq / 2).")
            return

        if (numTotalSample < 0 or numShftSample < 0):
            self.ShowErr("Total number of samples and sample shift number both need to be integers, greater than zero.")
            return

        self.__DesignFilter()

    def _varWindowTypeChanged(self, event):
 
        self.__stWinTyp = self._varWindowType.get()
        pv(self.__stWinTyp)

        self.__ComputeTimeVec()
        self.__ComputeFreqVec()
        self.__ComputeWindow()
        self.__ComputeWindowDFT()

        self._WindowFunctioninTimeDomain.update_yData(0, self.__window)
        self._WindowFunctioninTimeDomain.draw()

        self._WindowFunctionFrequencyResponse.update_yData(0, self.__winMag)
        self._WindowFunctionFrequencyResponse.draw()		

    def _info(self):
        self.show_info("Info", "FIR Filter Designer\nWritten by Philip M. Salmony\n29 November 2019\nphilsal.co.uk")

    def _varFilterTypeChanged(self):
        stFiltTyp = FilrType(self._varFilterType.get())
        pv(stFiltTyp)

        if (stFiltTyp == FilrType.LowPass):
            self._txtLowFrequency.configure(state='normal')
            self._txtHighFrequency.configure(state='disabled')
        elif (stFiltTyp == FilrType.HighPass):
            self._txtLowFrequency.configure(state='disabled')
            self._txtHighFrequency.configure(state='normal')
        elif (stFiltTyp == FilrType.BandPass):
            self._txtLowFrequency.configure(state='normal')
            self._txtHighFrequency.configure(state='normal')
        elif (stFiltTyp == FilrType.BandStop):
            self._txtLowFrequency.configure(state='normal')
            self._txtHighFrequency.configure(state='normal')

    def _export_coefficients(self):

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

    def _export_time_domain_data(self):

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

    def _export_frequency_domain_data(self):

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
    proj_path = os.path.dirname(os.path.abspath(__file__))
    if getattr(sys, 'frozen', False):
        # po("script is packaged!")
        proj_path = os.path.dirname(os.path.abspath(sys.executable))
    proj_path = os.path.join(proj_path, "..")

    app = FIRDesignerApp()
    win_xml = os.path.join(proj_path, "resources", 'window.xml')
    app.create_window(win_xml)
    menu_xml = os.path.join(proj_path, "resources", 'menu.xml')
    app.config_menu(menu_xml)
    app.init_paras()
    app.go()
