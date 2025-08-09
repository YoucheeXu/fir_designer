#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import sys
import os
import math
from enum import IntEnum, unique
from typing import override, cast

import numpy as np
from tkinter import filedialog as tkFileDialog
# import tkinter as tk

from pyutilities.logit import pv
from pyutilities.tkwin import EntryCtrl, RadioButtonGroupCtrl, ComboboxCtrl
from pyutilities.tkwin import tkWin
from pyutilities.matplot import MatPlotCtrl
from pyutilities.matplot import LineData


_ = np.seterr(divide='ignore', invalid='ignore')


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
class FilrType(IntEnum):
    LowPass = 0
    HighPass = 1
    BandPass = 2
    BandStop = 3


@unique
class RespPltType(IntEnum):
    Impulse = 0
    Step = 1


class FIRDesignerApp(tkWin):

    def __init__(self, cur_path: str, xmlfile: str):
        super().__init__(cur_path, xmlfile)
        self._win_typ: str = ""

        # Initial settings
        self._freq_sample_num: int = 256

        self._txt_sample_time: EntryCtrl = cast(EntryCtrl, self.get_control("SampleTime"))
        self._txt_sample_time.set_val("0.01")
        self._txt_low_freq: EntryCtrl = cast(EntryCtrl, self.get_control("LowFrequency"))
        self._txt_low_freq.set_val("10.0")
        self._txt_high_freq: EntryCtrl = cast(EntryCtrl, self.get_control("HighFrequency"))
        self._txt_high_freq.set_val("20.0")
        self._txt_filter_len: EntryCtrl = cast(EntryCtrl, self.get_control("FilterLength"))
        self._txt_filter_len.set_val("64")
        self._txt_shift_samples: EntryCtrl = cast(EntryCtrl, self.get_control("ShiftSamples"))
        self._txt_shift_samples.set_val("32")
        self._cmb_win_typ: ComboboxCtrl = cast(ComboboxCtrl, self.get_control("WindowType"))
        self._cmb_win_typ.set_val("Hamming")
        self._rad_filt_typ: RadioButtonGroupCtrl = cast(RadioButtonGroupCtrl, self.get_control("FilterType"))
        self._rad_filt_typ.set_val(0)    # FilrType.LowPass
        self._rad_resp_plot: RadioButtonGroupCtrl = cast(RadioButtonGroupCtrl, self.get_control("ResponsePlot"))
        self._rad_resp_plot.set_val(0)    # Impulse
        self._cmb_display_typ: ComboboxCtrl = cast(ComboboxCtrl, self.get_control("DisplayType"))
        self._cmb_display_typ.set_val("Both")        

        # Time Domain
        total_sample_num = int(self._txt_filter_len.get_val())
        self._time_vec: list[float] = [0.0] * total_sample_num
        self._impulse_resp: list[float] = [0.0] * total_sample_num
        self._step_resp: list[float] = [0.0] * total_sample_num
        self._window: list[float] = [0.0] * total_sample_num
        self._wind_imp_resp: list[float] = [0.0] * total_sample_num
        self._wind_step_resp: list[float] = [0.0] * total_sample_num

        self._filter_resp_time_domain: MatPlotCtrl = \
            cast(MatPlotCtrl, self.get_control("FilterResponseinTimeDomain"))
        self._win_funct_time_domain: MatPlotCtrl = \
            cast(MatPlotCtrl, self.get_control("WindowFunctioninTimeDomain"))

        self._filter_freq_resp: MatPlotCtrl = \
            cast(MatPlotCtrl, self.get_control("FilterFrequencyResponse"))
        self._win_funct_freq_resp: MatPlotCtrl = \
            cast(MatPlotCtrl, self.get_control("WindowFunctionFrequencyResponse"))

        # Frequency Domain
        freq_sample_num = self._freq_sample_num
        self._freq_vec: list[float] = [0.0] * freq_sample_num
        self._imp_resp_mag: list[float] = [0.0] * freq_sample_num
        self._win_resp_mag: list[float] = [0.0] * freq_sample_num
        self._win_mag: list[float] = [0.0] * freq_sample_num

    def _before_go(self):
        self._rad_filt_type_changed()
        self._design_filter()

    def _design_filter(self):
        self._compute_time_vec()
        self._compute_win()
        self._compute_resps()
        self._compute_wind_resps()

        self._compute_freq_vec()
        self._compute_resp_bode()
        self._compute_win_dft()

        self._update_charts()
        self._update_plot_settings()

    def _update_charts(self):

        self._filter_resp_time_domain.xdata = self._time_vec
        y_line = LineData(self._impulse_resp, {"label": "Impulse Response"})
        _ = self._filter_resp_time_domain.add_line(y_line)
        y_line = LineData(self._wind_imp_resp, {"label": "Windowed Impulse Response"})
        _ = self._filter_resp_time_domain.add_line(y_line)
        y_line = LineData(self._step_resp, {"label": "Step Response"}, visible = False)
        _ = self._filter_resp_time_domain.add_line(y_line)
        y_line = LineData(self._wind_step_resp, {"label": "Windowed Step Response"}, visible = False)
        _ = self._filter_resp_time_domain.add_line(y_line)
        self._filter_resp_time_domain.draw()

        self._win_funct_time_domain.xdata = self._time_vec
        y_line = LineData(self._window)
        _ = self._win_funct_time_domain.add_line(y_line)
        self._win_funct_time_domain.draw()

        self._filter_freq_resp.xdata = self._freq_vec
        y_line = LineData(self._imp_resp_mag)
        id_ = self._filter_freq_resp.add_line(y_line)
        print(f"id of imp_resp_mag: {id_}")
        print(f"length of _imp_resp_mag: {len(self._imp_resp_mag)}")
        y_line = LineData(self._win_resp_mag)
        id_ = self._filter_freq_resp.add_line(y_line)
        print(f"id of _win_resp_mag: {id_}")
        print(f"length of _win_resp_mag: {len(self._win_resp_mag)}")
        self._filter_freq_resp.draw()

        self._win_funct_freq_resp.xdata = self._freq_vec
        y_line = LineData(self._win_mag)
        _ = self._win_funct_freq_resp.add_line(y_line)
        self._win_funct_freq_resp.draw()

    # Time Domain Functions
    def _compute_time_vec(self):
        total_sample_num = int(self._txt_filter_len.get_val())
        self._time_vec = [0.0] * total_sample_num

        ti_sample = float(self._txt_sample_time.get_val())
        for n in range(total_sample_num):
            self._time_vec[n] = n * ti_sample

        # pv(self._timeVec)

    def _compute_resps(self):
        total_sample_num = int(self._txt_filter_len.get_val())
        self._impulse_resp = [0.0] * total_sample_num
        self._step_resp = [0.0] * total_sample_num

        shft_sample_num = int(self._txt_shift_samples.get_val())
        filt_typ = FilrType(self._rad_filt_typ.get_val())
        ti_sample = float(self._txt_sample_time.get_val())
        frq_cutoff_lo = float(self._txt_low_freq.get_val())
        frq_cutoff_hi = float(self._txt_high_freq.get_val())
        for n in range(total_sample_num):
            if n != shft_sample_num:
                if filt_typ == FilrType.LowPass:
                    self._impulse_resp[n] = math.sin(2.0 * math.pi * frq_cutoff_lo * \
                        ti_sample * (n - shft_sample_num)) / (math.pi * ti_sample * (n - shft_sample_num))
                elif filt_typ == FilrType.HighPass:
                    self._impulse_resp[n] = (math.sin(math.pi * (n - shft_sample_num)) - \
                        math.sin(2.0 * math.pi * frq_cutoff_hi * ti_sample * (n - shft_sample_num))) \
                        / (math.pi * ti_sample * (n - shft_sample_num))
                elif filt_typ == FilrType.BandPass:
                    self._impulse_resp[n] = (math.sin(2.0 * math.pi * frq_cutoff_hi * \
                        ti_sample * (n - shft_sample_num)) - math.sin(2.0 * math.pi * \
                        frq_cutoff_lo * ti_sample * (n - shft_sample_num))) / (math.pi * \
                        ti_sample * (n - shft_sample_num))
                elif filt_typ == FilrType.BandStop:
                    self._impulse_resp[n] = (math.sin(2.0 * math.pi * frq_cutoff_lo * \
                        ti_sample * (n - shft_sample_num)) - math.sin(2.0 * math.pi * \
                        frq_cutoff_hi * ti_sample * (n - shft_sample_num)) + math.sin(math.pi * \
                        (n - shft_sample_num))) / (math.pi * ti_sample * (n - shft_sample_num))

            else: # Avoid divide-by-zero, limit is 2*fc
                if filt_typ == FilrType.LowPass:
                    self._impulse_resp[n] = 2.0 * frq_cutoff_lo
                elif filt_typ == FilrType.HighPass:
                    self._impulse_resp[n] = 1.0 / ti_sample - 2.0 * frq_cutoff_hi
                elif filt_typ == FilrType.BandPass:
                    self._impulse_resp[n] = 2.0 * frq_cutoff_hi - 2.0 * frq_cutoff_lo
                elif filt_typ == FilrType.BandStop:
                    self._impulse_resp[n] = 2.0 * frq_cutoff_lo - 2.0 * frq_cutoff_hi + 1.0 / ti_sample

        # Normalise by DC gain to achieve 0dB gain at DC and then compute step response
        for n in range(total_sample_num):
            self._impulse_resp[n] *= ti_sample

            if n == 0:
                self._step_resp[n] = 0.5 * self._impulse_resp[n]
            else:
                self._step_resp[n] = self._step_resp[n - 1] + 0.5 * (self._impulse_resp[n] + \
                    self._impulse_resp[n - 1])

        # pv(self._impulse_resp)
        # pv(self._step_resp)

    def _compute_win(self):
        total_sample_num = int(self._txt_filter_len.get_val())
        self._window = [0.0] * total_sample_num

        win_typ = self._cmb_win_typ.get_val()

        for n in range(total_sample_num):
            if win_typ == "Rectangular":
                self._window[n] = 1.0
            elif win_typ == "Triangular":
                self._window[n] = 1.0 - abs((n - 0.5 * total_sample_num) / (0.5 * total_sample_num))
            elif win_typ == "Welch":
                self._window[n] = 1.0 - math.pow((n - 0.5 * total_sample_num) / (0.5 * total_sample_num), 2.0)
            elif win_typ == "Sine":
                self._window[n] = math.sin(math.pi * n / (total_sample_num))
            elif win_typ == "Hann":
                self._window[n] = 0.5 * (1 - math.cos(2.0 * math.pi * n / (total_sample_num)))
            elif win_typ == "Hamming":
                self._window[n] = (25.0 / 46.0) - (21.0 / 46.0) * \
                    math.cos(2.0 * math.pi * n / (total_sample_num))
            elif win_typ == "Blackman":
                self._window[n] = 0.42 - 0.5 * math.cos(2.0 * math.pi * \
                    n / (total_sample_num)) + 0.08 * math.cos(4.0 * math.pi * n / (total_sample_num))
            elif win_typ == "Nuttall":
                self._window[n] = 0.355768 - 0.487396 * math.cos(2.0 * math.pi * \
                    n / (total_sample_num)) + 0.144232 * math.cos(4.0 * math.pi * \
                    n / (total_sample_num)) - 0.012604 * math.cos(6.0 * math.pi * n / (total_sample_num))
            elif win_typ == "BlackmanNuttall":
                self._window[n] = 0.3635819 - 0.4891775 * math.cos(2.0 * math.pi * \
                    n / (total_sample_num)) + 0.1365995 * math.cos(4.0 * math.pi * \
                    n / (total_sample_num)) - 0.0106411 * math.cos(6.0 * math.pi * n / (total_sample_num))
            elif win_typ == "BlackmanHarris":
                self._window[n] = 0.35875 - 0.48829 * math.cos(2.0 * math.pi * \
                    n / (total_sample_num)) + 0.14128 * math.cos(4.0 * math.pi * \
                    n / (total_sample_num)) - 0.01168 * math.cos(6.0 * math.pi * n / (total_sample_num))
            elif win_typ == "FlatTop":
                self._window[n] = 0.21557895 - 0.41663158 * math.cos(2.0 * math.pi * \
                    n / (total_sample_num)) + 0.277263158 * math.cos(4.0 * math.pi * \
                    n / (total_sample_num)) - 0.083578947 * math.cos(6.0 * math.pi * \
                    n / (total_sample_num)) + 0.006947368 * math.cos(8.0 * math.pi * n / (total_sample_num))
            else:
                self._window[n] = 1.0

        # pv(self._window)

    def _compute_wind_resps(self):
        total_sample_num = int(self._txt_filter_len.get_val())
        self._wind_imp_resp = [0.0] * total_sample_num
        self._wind_step_resp = [0.0] * total_sample_num

        for n in range(total_sample_num):
            self._wind_imp_resp[n] = self._impulse_resp[n] * self._window[n]

            if n == 0:
                self._wind_step_resp[n] = 0.5 * self._wind_step_resp[n]
            else:
                self._wind_step_resp[n] = self._wind_step_resp[n - 1] + \
                    0.5 * (self._wind_imp_resp[n] + self._wind_imp_resp[n - 1])

        # pv(self._wind_imp_resp)
        # pv(self._wind_step_resp)

    def _compute_freq_vec(self):
        ''' Frequency Domain Functions
        '''
        freq_sample_num = self._freq_sample_num
        self._freq_vec = [0.0] * freq_sample_num

        ti_sample = float(self._txt_sample_time.get_val())
        df = (0.5 / ti_sample) / (freq_sample_num - 1.0)
        pv(df)

        for n in range(freq_sample_num):
            self._freq_vec[n] = n * df

    def _compute_resp_bode(self):
        freq_sample_num = self._freq_sample_num
        self._imp_resp_mag = [0.0] * freq_sample_num
        self._win_resp_mag = [0.0] * freq_sample_num

        numTotalSample = int(self._txt_filter_len.get_val())
        ti_sample = float(self._txt_sample_time.get_val())
        for fIndex in range(freq_sample_num):
            re = 0.0
            im = 0.0
            re_win = 0.0
            im_win = 0.0

            for n in range(numTotalSample):
                re = re + self._impulse_resp[n] * math.cos(2.0 * math.pi * \
                    self._freq_vec[fIndex] * n * ti_sample)
                im = im - self._impulse_resp[n] * math.sin(2.0 * math.pi * \
                    self._freq_vec[fIndex] * n * ti_sample)
                re_win = re_win + self._wind_imp_resp[n] * math.cos(2.0 * math.pi * \
                    self._freq_vec[fIndex] * n * ti_sample)
                im_win = im_win - self._wind_imp_resp[n] * math.sin(2.0 * math.pi * \
                    self._freq_vec[fIndex] * n * ti_sample)

            self._imp_resp_mag[fIndex] = 10.0 * math.log10(re * re + im * im)
            self._win_resp_mag[fIndex] = 10.0 * math.log10(re_win * re_win + im_win * im_win)

    def _get_gain_cutoff(self) -> float:
        re = 0.0
        im = 0.0

        total_sample_num = int(self._txt_filter_len.get_val())
        ti_sample = float(self._txt_sample_time.get_val())
        frq_cutoff_lo = float(self._txt_low_freq.get_val())
        # frq_cutoff_hi = float(self._txt_high_freq.get_val())
        for n in range(total_sample_num):
            re = re + self._impulse_resp[n] * math.cos(2.0 * math.pi * \
                frq_cutoff_lo * n * ti_sample)
            im = im - self._impulse_resp[n] * math.sin(2.0 * math.pi * \
                frq_cutoff_lo * n * ti_sample)

        return (10.0 * math.log10(re * re + im * im))

    def _compute_win_dft(self):
        freq_sample_num = self._freq_sample_num
        self._win_mag = [0.0] * freq_sample_num

        total_sample_num = int(self._txt_filter_len.get_val())
        ti_sample = float(self._txt_sample_time.get_val())
        for i in range(freq_sample_num):
            re = 0.0
            im = 0.0
            for n in range(total_sample_num):
                re = re + self._window[n] * math.cos(2.0 * math.pi * \
                    self._freq_vec[i] * n * ti_sample)
                im = im - self._window[n] * math.sin(2.0 * math.pi * \
                    self._freq_vec[i] * n * ti_sample)

            self._win_mag[i] = 10.0 * math.log10(re * re + im * im)

    def _update_plot_settings(self):
        resp_plt_typ = RespPltType(self._rad_resp_plot.get_val())
        pv(resp_plt_typ.name)
        display_typ = self._cmb_display_typ.get_val()
        pv(display_typ)

        if resp_plt_typ == RespPltType.Impulse:        # Impulse
            if display_typ == "Both":
                self._filter_resp_time_domain.show_line(0, True)
                self._filter_resp_time_domain.show_line(1, True)
            elif display_typ == "Raw Only":
                self._filter_resp_time_domain.show_line(0, True)
                self._filter_resp_time_domain.show_line(1, False)
            elif display_typ == "Win Only":
                self._filter_resp_time_domain.show_line(0, False)
                self._filter_resp_time_domain.show_line(1, True)                
            self._filter_resp_time_domain.show_line(2, False)
            self._filter_resp_time_domain.show_line(3, False)            

        else:        # Step
            if display_typ == "Both":
                self._filter_resp_time_domain.show_line(2, True)
                self._filter_resp_time_domain.show_line(3, True)
            elif display_typ == "Raw Only":
                self._filter_resp_time_domain.show_line(2, True)
                self._filter_resp_time_domain.show_line(3, False)
            elif display_typ == "Win Only":
                self._filter_resp_time_domain.show_line(2, False)
                self._filter_resp_time_domain.show_line(3, True)
            self._filter_resp_time_domain.show_line(0, False)
            self._filter_resp_time_domain.show_line(1, False)

        if display_typ == "Both":
            self._filter_freq_resp.show_line(0, True)
            self._filter_freq_resp.show_line(1, True)
        elif display_typ == "Raw Only":
            self._filter_freq_resp.show_line(0, True)
            self._filter_freq_resp.show_line(1, False)
        elif display_typ == "Win Only":
            self._filter_freq_resp.show_line(0, False)
            self._filter_freq_resp.show_line(1, True)

        self._filter_resp_time_domain.draw()
        self._filter_freq_resp.draw()

    def _btn_designfilter_click(self):
        ti_sample = float(self._txt_sample_time.get_val())
        pv(ti_sample)
        total_sample_num = int(self._txt_filter_len.get_val())
        pv(total_sample_num)
        shft_sample_num = int(self._txt_shift_samples.get_val())
        pv(shft_sample_num)
        win_typ = self._cmb_win_typ.get_val()
        pv(win_typ)
        filt_typ = FilrType(self._rad_filt_typ.get_val())
        pv(filt_typ)

        frq_cutoff_lo = 0
        frq_cutoff_hi = 0
        if filt_typ == FilrType.LowPass:
            frq_cutoff_lo = float(self._txt_low_freq.get_val())
            pv(frq_cutoff_lo)
        elif filt_typ == FilrType.HighPass:
            frq_cutoff_hi = float(self._txt_high_freq.get_val())
            pv(frq_cutoff_hi)
        elif filt_typ == FilrType.BandPass:
            frq_cutoff_lo = float(self._txt_low_freq.get_val())
            frq_cutoff_hi = float(self._txt_high_freq.get_val())
            pv(frq_cutoff_lo)
            pv(frq_cutoff_hi)
        elif filt_typ == FilrType.BandStop:
            frq_cutoff_lo = float(self._txt_low_freq.get_val())
            frq_cutoff_hi = float(self._txt_high_freq.get_val())
            pv(frq_cutoff_lo)
            pv(frq_cutoff_hi)

        if ti_sample < 0.0:
            self.show_err("Sampling frequency cannot be negative.")
            return

        if frq_cutoff_lo >= 0.5 / ti_sample or frq_cutoff_hi >= 0.5 / ti_sample:
            self.show_err("Cut-off frequency has to be less than the " + \
                "Nyquist frequency (i.e. sampling freq / 2).")
            return

        if total_sample_num < 0 or shft_sample_num < 0:
            self.show_err("Total number of samples and sample shift " + \
                "number both need to be integers, greater than zero.")
            return

        self._design_filter()

    def _var_win_typ_changed(self):
 
        self._win_typ = self._cmb_win_typ.get_val()
        pv(self._win_typ)

        self._compute_time_vec()
        self._compute_freq_vec()
        self._compute_win()
        self._compute_win_dft()

        self._win_funct_time_domain.update_ydata(0, self._window)
        self._win_funct_time_domain.draw()

        self._win_funct_freq_resp.update_ydata(0, self._win_mag)
        self._win_funct_freq_resp.draw()        

    def _info(self):
        self.show_info("Info", "FIR Filter Designer\nWritten " + \
            "by Philip M. Salmony\n29 November 2019\nphilsal.co.uk")

    def _rad_filt_type_changed(self):
        filt_typ = FilrType(self._rad_filt_typ.get_val())
        pv(filt_typ)

        if filt_typ == FilrType.LowPass:
            _ = self._txt_low_freq.configure(state='normal')
            _ = self._txt_high_freq.configure(state='disabled')
        elif filt_typ == FilrType.HighPass:
            _ = self._txt_low_freq.configure(state='disabled')
            _ = self._txt_high_freq.configure(state='normal')
        elif filt_typ == FilrType.BandPass:
            _ = self._txt_low_freq.configure(state='normal')
            _ = self._txt_high_freq.configure(state='normal')
        elif filt_typ == FilrType.BandStop:
            _ = self._txt_low_freq.configure(state='normal')
            _ = self._txt_high_freq.configure(state='normal')

    @override
    def process_message(self, idmsg: str, **kwargs):
        match idmsg:
            case "ResponsePlot":
                self._update_plot_settings()
            case "DisplayType":
                self._update_plot_settings()
            case "WindowType":
                self._var_win_typ_changed()
            case "FilterType":
                self._rad_filt_type_changed()
            case "DesignFilter":
                self._btn_designfilter_click()
            case "Info":
                self._info()
            case "ExportCoefficients":
                self._export_coefficients()
            case "ExportTimeDomainData":
                self._export_time_domain_data()
            case "ExportFrequencyDomainData":
                self._export_frequency_domain_data()
            case _:
                _ = super().process_message(idmsg, **kwargs)

    def _export_coefficients(self):

        file_types = [('Text File', '*.txt'), ('All File', '*.*')]
        title = "Export Filter Coefficients"

        r = tkFileDialog.asksaveasfilename(title=title, filetypes=file_types)

        if r != "":
            pv(r)
            ext = os.path.splitext(r)[-1]
            if not ext:
                r += ".txt"

            with open(r, 'w') as file:
                total_sample_num = int(self._txt_filter_len.get_val())
                ti_sample = float(self._txt_sample_time.get_val())
                frq_cutoff_lo = float(self._txt_low_freq.get_val())
                frq_cutoff_hi = float(self._txt_high_freq.get_val())

                # string[] data = new string[3]
                data_list = [""] * 3
                # dataLst[0] = "Filr Order: " + self._numTotalSample + \
                    # " Sampling Freq (Hz): " + (1.0 / self._tiSample).ToString("F6") + \
                    # " Cut-Off Freq Lo (Hz): " + self._frqCutOffLo.ToString("F6") + \
                    # " Cut-Off Freq Hi (Hz): " + self._frqCutOffHi.ToString("F6") + "\n\n"
                data_list[0] = (f"Filter Order: {total_sample_num} "\
                    f"Frequency: {1.0 / ti_sample}Hz, Cut-Off Frequency Low: " \
                    f"{frq_cutoff_lo}Hz, Cut-Off Frequency High: {frq_cutoff_hi}Hz\n\n")

                data_list[1] = f"Coefficient, {self._wind_imp_resp[0]}"    # F7
                data_list[2] = "float Coefficient[] = {" + str(self._wind_imp_resp[0]) + "f"    # F7

                for n in range(total_sample_num):
                    data_list[1] += "," + str(self._wind_imp_resp[n])            # F9
                    data_list[2] += "," + str(self._wind_imp_resp[n]) + "f"    # F7

                # dataLst[1] = ", ".join()
                # dataLst[2]

                data_list[1] += "\n\n"
                data_list[2] += "}"

                for data in data_list:
                    _ = file.write(data)

                self.show_info(title, f"Coefficients written to {r}")

    def _export_time_domain_data(self):

        filetypes = [('Text File', '*.txt'), ('All File', '*.*')]
        title = "Export Time Domain Data"

        if r := tkFileDialog.asksaveasfilename(title=title, filetypes=filetypes):
            pv(r)
            ext = os.path.splitext(r)[-1]
            if not ext:
                r += ".txt"

            with open(r, 'w') as file:
                total_sample_num = int(self._txt_filter_len.get_val())
                ti_sample = float(self._txt_sample_time.get_val())
                frq_cutoff_lo = float(self._txt_low_freq.get_val())
                frq_cutoff_hi = float(self._txt_high_freq.get_val())

                data_list = [""] * 4
                # dataLst[0] = "[TIME DOMAIN DATA (TIME/IMPULSE/STEP)] Filter Order: " + \
                    # "self._numTotalSample + " Sampling Frequency (Hz): " + \
                    # (1.0 / self._tiSample).ToString("F6") + \
                    # " Cut-Off Frequency Low (Hz): " + self._frqCutOffLo.ToString("F6") + \
                    # " Cut-Off Frequency High (Hz): " + self._frqCutOffHi.ToString("F6") + "\n\n"
                data_list[0] = (f"[TIME DOMAIN DATA (TIME/IMPULSE/STEP)]\n\nFilter Order: " \
                    f"{total_sample_num}, Sampling Frequency: {1.0 / ti_sample}Hz, " \
                    f"Cut-Off Frequency Low: {frq_cutoff_lo}Hz, Cut-Off Frequency High: {frq_cutoff_hi}Hz\n\n")

                data_list[1] = f"Time Vector (s), {self._time_vec[0]}"                  # F6
                data_list[2] = f"Windowed Impulse Response, {self._wind_imp_resp[0]}"   # F9
                data_list[3] = f"Windowed Step Response, {self._wind_step_resp[0]}"     # F9
                for n in range(total_sample_num):
                    data_list[1] += "," + str(self._time_vec[n])            # F6
                    data_list[2] += "," + str(self._wind_imp_resp[n])       # F9
                    data_list[3] += "," + str(self._wind_step_resp[n])      # F9
                # dataLst[1] = ", ".join(str(self._timeVec))
                # dataLst[2] = ", ".join(str(self._windowedImpResp))
                # dataLst[3] = ", ".join(str(self._windowedStepResp))

                data_list[1] += "\n\n"
                data_list[2] += "\n\n"

                for data in data_list:
                    _ = file.write(data)
                self.show_info(title, f"Data written to {r}")

    def _export_frequency_domain_data(self):

        filetypes = [('Text File', '*.txt'), ('All File', '*.*')]
        title = "Export Frequency Domain Data"

        if r := tkFileDialog.asksaveasfilename(title=title, filetypes=filetypes):
            pv(r)
            ext = os.path.splitext(r)[-1]
            if not ext:
                r += ".txt"

            with open(r, 'w') as file:
                total_sample_num = int(self._txt_filter_len.get_val())
                ti_sample = float(self._txt_sample_time.get_val())
                frq_cutoff_lo = float(self._txt_low_freq.get_val())
                frq_cutoff_hi = float(self._txt_high_freq.get_val())

                data_list = [""] * 4
                # dataLst[0] = "[FREQUENCY DOMAIN DATA (FREQUENCY/RAW/WINDOWED)] Filr Order: " + \
                    # self._numTotalSample + " Sampling Freq (Hz): " + \
                    # (1.0 / self._tiSample).ToString("F6") + \
                    # " Cut-Off Freq Lo (Hz): " + self._frqCutOffLo.ToString("F6") + \
                    # " Cut-Off Freq Hi (Hz): " + self._frqCutOffHi.ToString("F6") + "\n\n"
                data_list[0] = (f"[FREQUENCY DOMAIN DATA (FREQUENCY/RAW/WINDOWED)]\n\n"
                    f"Filter Order: {total_sample_num}, Sampling Frequency: {1.0 / ti_sample}"
                    f"Hz, Cut-Off Frequency Low: {frq_cutoff_lo}Hz, Cut-Off Frequency High: "
                    f"{frq_cutoff_hi}Hz\n\n")

                data_list[1] = f"Frequency Vecctor (Hz), {self._freq_vec[0]}"               # F6
                data_list[2] = f"Impulse Response Magnitude, {self._imp_resp_mag[0]}"       # F9
                data_list[3] = f"Window Response Magnitude, {self._win_resp_mag[0]}"        # F9
                for n in range(self._freq_sample_num):
                    data_list[1] += "," + str(self._freq_vec[n])        # F6
                    data_list[2] += "," + str(self._imp_resp_mag[n])    # F9
                    data_list[3] += "," + str(self._win_resp_mag[n])    # F9

                data_list[1] += "\n\n"
                data_list[2] += "\n\n"

                for data in data_list:
                    _ = file.write(data)
                self.show_info(title, f"Data written to {r}")


if __name__ == '__main__':
    cur_path = os.path.dirname(os.path.abspath(__file__))
    if getattr(sys, 'frozen', False):
        # print("script is packaged!")
        cur_path = os.path.dirname(os.path.abspath(sys.executable))
    proj_path = os.path.join(cur_path, "..")
    win_xml = os.path.join(proj_path, 'resources', 'window.xml')
    my_app = FIRDesignerApp(proj_path, win_xml)
    my_app.go()
