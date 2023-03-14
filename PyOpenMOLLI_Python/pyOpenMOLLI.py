# Import functions
import copy
from math import pi, sqrt, ceil, floor
from types import SimpleNamespace
import types
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts
from pypulseq.make_extended_trapezoid import make_extended_trapezoid
from pypulseq.make_trig import make_trig as make_trigger  # For trigger


def bssfp_readout(seq, system, fov=200e-3, Nstartup=11, Ny=128):
    """
    Creates a Balanced steady-state free precession (bSSFP) sequence and adds
    to seq object

    Parameters
    ----------
    seq: object
        Sequence object
    system : Opts
        System limits.
    fov : float, optional
        Field-of-view [m]. Default is 0.2 m.
    Nstartup : float, optional
        Number of start-up RF pulses. Default is 11
    Ny : float, optional
        Number of phase encoding lines. Default is 128.
   Returns
    -------
    seq : SimpleNamespace
        Seq object with bSSFP readout
    TR : float
        Repetition time
    Ny : float
        Final number of phase encoding lines.
    """
    # Sequence Parameters
    enc = 'xyz'
    Nx = 128
    thk = 6e-3
    fa = 35  # [deg]

    Nramp = Nstartup

    # ADC duration (controls TR/TE)
    adc_dur = 2560 / 2  # [us]

    rf_dur = 490  # [us]
    rf_apo = 0.5
    rf_bwt = 1.5

    #############################################################################
    #                 Create slice selection pulse and gradient
    rf, g_ss, __ = make_sinc_pulse(flip_angle=fa * pi / 180, system=system, duration=rf_dur * 1e-6, slice_thickness=thk,
                                   apodization=rf_apo, time_bw_product=rf_bwt)  # for next pyPulseq version add:, return_gz=True)
    g_ss.channel = enc[2]

    # Slice refocusing
    g_ss_reph = make_trapezoid(channel=enc[2], system=system, area=-g_ss.area / 2, duration=0.00017 * 2)
    #############################################################################

    rf.delay = calc_duration(g_ss) - calc_duration(rf) + rf.delay

    #############################################################################
    #                        Readout gradient and ADC
    delta_k = 1 / fov
    kWidth = Nx * delta_k

    # Readout and ADC
    g_ro = make_trapezoid(channel=enc[0], system=system, flat_area=kWidth, flat_time=(adc_dur) * 1e-6)
    adc = make_adc(num_samples=Nx, duration=g_ro.flat_time, delay=g_ro.rise_time)

    # Readout rewinder
    g_ro_pre = make_trapezoid(channel=enc[0], system=system, area=-g_ro.area / 2)

    phaseAreas_tmp = np.arange(0, Ny, 1).tolist()
    phaseAreas = np.dot(np.subtract(phaseAreas_tmp, Ny / 2), delta_k)

    gs8_times = [0, g_ro.fall_time,
                 g_ro.fall_time + g_ro_pre.rise_time,
                 g_ro.fall_time + g_ro_pre.rise_time + g_ro_pre.flat_time,
                 g_ro.fall_time + g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time]
    gs8_amp = [g_ro.amplitude, 0, g_ro_pre.amplitude, g_ro_pre.amplitude, 0]
    gx_2 = make_extended_trapezoid(channel=enc[0], times=gs8_times, amplitudes=gs8_amp)

    # Calculate phase encoding gradient duration
    pe_dur = calc_duration(gx_2)

    gx_allExt_times = [0, g_ro_pre.rise_time, g_ro_pre.rise_time + g_ro_pre.flat_time,
                       g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time,
                       g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time + g_ro.rise_time,
                       g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time + g_ro.rise_time + g_ro.flat_time + 1e-5,
                       g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time + g_ro.rise_time + g_ro.flat_time + 1e-5
                       + g_ro.fall_time,
                       g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time + g_ro.rise_time + g_ro.flat_time + 1e-5
                       + g_ro.fall_time + g_ro_pre.rise_time,
                       g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time + g_ro.rise_time + g_ro.flat_time + 1e-5
                       + g_ro.fall_time + g_ro_pre.rise_time + g_ro_pre.flat_time,
                       g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time + g_ro.rise_time + g_ro.flat_time + 1e-5
                       + g_ro.fall_time + g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time]
    gx_allExt_amp = [0, g_ro_pre.amplitude, g_ro_pre.amplitude, 0, g_ro.amplitude, g_ro.amplitude, 0,
                     g_ro_pre.amplitude, g_ro_pre.amplitude, 0]
    gx_all = make_extended_trapezoid(channel=enc[0], times=gx_allExt_times, amplitudes=gx_allExt_amp)
    #############################################################################

    gzrep_times = [0, g_ss_reph.rise_time,
                   g_ss_reph.rise_time + g_ss_reph.flat_time,
                   g_ss_reph.rise_time + g_ss_reph.flat_time + g_ss_reph.fall_time,
                   g_ss_reph.rise_time + g_ss_reph.flat_time + g_ss_reph.fall_time + calc_duration(g_ro)
                   + 2 * calc_duration(g_ro_pre) - 2 * calc_duration(g_ss_reph) + 1e-5,
                   g_ss_reph.rise_time + g_ss_reph.flat_time + g_ss_reph.fall_time + calc_duration(g_ro)
                   + 2 * calc_duration(g_ro_pre) - 2 * calc_duration(g_ss_reph) + 1e-5 + g_ss_reph.rise_time,
                   g_ss_reph.rise_time + g_ss_reph.flat_time + g_ss_reph.fall_time + calc_duration(g_ro)
                   + 2 * calc_duration(g_ro_pre) - 2 * calc_duration(
                       g_ss_reph) + 1e-5 + g_ss_reph.rise_time + g_ss_reph.flat_time,
                   g_ss_reph.rise_time + g_ss_reph.flat_time + g_ss_reph.fall_time + calc_duration(g_ro)
                   + 2 * calc_duration(g_ro_pre) - 2 * calc_duration(
                       g_ss_reph) + 1e-5 + g_ss_reph.rise_time + g_ss_reph.flat_time + g_ss_reph.fall_time]
    gzrep_amp = [0, g_ss_reph.amplitude, g_ss_reph.amplitude, 0, 0, g_ss_reph.amplitude, g_ss_reph.amplitude, 0]
    gzrep_all = make_extended_trapezoid(channel=enc[2], times=gzrep_times, amplitudes=gzrep_amp)

    adc.delay = g_ro_pre.rise_time + g_ro_pre.flat_time + g_ro_pre.fall_time + g_ro.rise_time + 0.5e-5

    # finish timing calculation
    TR = calc_duration(g_ss) + calc_duration(gx_all)
    TE = TR / 2

    ni_acqu_pattern = np.arange(1, Ny + 1, 1)
    Ny_aq = len(ni_acqu_pattern)
    print('Acquisition window is: %3.2f ms' % (TR * Ny_aq * 1e3))

    rf05 = rf
    rf_waveform = rf.signal
    ############################################################################
    #                          Start-up RF pulses
    #                          (ramp-up of Nramp)
    for nRamp in np.arange(1, Nramp + 1, 1):
        if np.mod(nRamp, 2):
            rf.phase_offset = 0
            adc.phase_offset = 0
        else:
            rf.phase_offset = -pi
            adc.phase_offset = -pi

        rf05.signal = np.divide(nRamp, Nramp) * rf_waveform

        gyPre_2 = make_trapezoid('y', area=phaseAreas[0], duration=pe_dur, system=system)
        gyPre_1 = make_trapezoid('y', area=-phaseAreas[0], duration=pe_dur, system=system)

        gyPre_times = [0, gyPre_2.rise_time,
                       gyPre_2.rise_time + gyPre_2.flat_time,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time + g_ro.flat_time + 1e-5,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time + g_ro.flat_time + 1e-5
                       + gyPre_1.rise_time,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time + g_ro.flat_time + 1e-5
                       + gyPre_1.rise_time + gyPre_1.flat_time,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time + g_ro.flat_time + 1e-5
                       + gyPre_1.rise_time + gyPre_1.flat_time + gyPre_1.fall_time]
        gyPre_amp = [0, gyPre_2.amplitude, gyPre_2.amplitude, 0, 0, gyPre_1.amplitude, gyPre_1.amplitude, 0]
        gyPre_all = make_extended_trapezoid(channel=enc[1], times=gyPre_times, amplitudes=gyPre_amp)

        if nRamp == 1:
            seq.add_block(rf05, g_ss)
            seq.add_block(gx_all, gyPre_all, gzrep_all)
        else:
            seq.add_block(rf05, g_ss)
            seq.add_block(gx_all, gyPre_all, gzrep_all)
    ############################################################################

    ############################################################################
    #           Actual Readout iterate number of phase encoding steps Ny
    for i in np.arange(1, Ny + 1, 1):
        # *******************************
        # Phase cycling
        if np.mod(i + Nramp, 2):
            rf.phase_offset = 0
            adc.phase_offset = 0
        else:
            rf.phase_offset = -pi
            adc.phase_offset = -pi
        # *******************************

        # ***************************************************************************************
        #                     Create phase encoding gradients
        if np.mod(Nramp, 2):
            gyPre_2 = make_trapezoid('y', area=phaseAreas[i - 1], duration=pe_dur, system=system)
            if i > 1:
                gyPre_1 = make_trapezoid('y', area=-phaseAreas[np.mod(i + Ny - 2, Ny)],
                                         duration=pe_dur, system=system)
            else:
                gyPre_1 = make_trapezoid('y', area=-phaseAreas[np.mod(i + Ny - 1, Ny)],
                                         duration=pe_dur, system=system)
        else:
            gyPre_2 = make_trapezoid('y', area=phaseAreas[i], duration=pe_dur, system=system)
            if i > 1:
                gyPre_1 = make_trapezoid('y', area=-phaseAreas[np.mod(i + Ny - 2, Ny)],
                                         duration=pe_dur, system=system)
            else:
                gyPre_1 = make_trapezoid('y', area=phaseAreas[np.mod(i + Ny - 2, Ny)],
                                         duration=pe_dur, system=system)

        gyPre_times = [0, gyPre_2.rise_time,
                       gyPre_2.rise_time + gyPre_2.flat_time,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time + g_ro.flat_time + 1e-5,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time + g_ro.flat_time + 1e-5
                       + gyPre_2.rise_time,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time + g_ro.flat_time + 1e-5
                       + gyPre_2.rise_time + gyPre_2.flat_time,
                       gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time + g_ro.flat_time + 1e-5
                       + gyPre_2.rise_time + gyPre_2.flat_time + gyPre_2.fall_time]
        gyPre_amp = [0, gyPre_2.amplitude, gyPre_2.amplitude, 0, 0, -gyPre_2.amplitude, -gyPre_2.amplitude, 0]

        # Verify if phase encoding gradient amplitude is not zero
        try:
            gyPre_all = make_extended_trapezoid(channel=enc[1], times=gyPre_times, amplitudes=gyPre_amp)
            flag_zerogy = 0
        except:
            gyPre_times = [0, pe_dur,
                           pe_dur + g_ro.flat_time + 1e-5,
                           pe_dur + g_ro.flat_time + 1e-5 + pe_dur]
            gyPre_amp = [0, gyPre_2.amplitude, 0, 0]
            flag_zerogy = 1
        # ***************************************************************************************

        # add RF and Slice selecitve gradient
        seq.add_block(rf, g_ss)

        # add Phase and frequency encoding gradients and ADC
        if flag_zerogy == 1:
            seq.add_block(gx_all, gzrep_all, adc)
        else:
            seq.add_block(gx_all, gyPre_all, gzrep_all, adc)

    return seq, TR, Ny


def make_hyperSec_pulse(system: Opts = Opts(), duration: float = 0, freq_offset: float = 0,
                        phase_offset: float = 0, center_pos: float = 0.5,
                        delay: float = 0):
    """
    Creates a radio-frequency hyperSec pulse event

    Parameters
    ----------
    system : Opts, optional
        System limits. Default is a system limits object initialised to default values.
    duration : float, optional
        Duration in milliseconds (ms). Default is 0.
    phase_offset : float, optional
        Phase offset in Hertz (Hz). Default is 0.
    center_pos : float, optional
        Position of peak. Default is 0.5
    delay : float, optional
        Delay in milliseconds (ms). Default is 0.
   Returns
    -------
    rf_inv : SimpleNamespace
        Radio-frequency inversion pulse event.
    """

    Beta = 674.1917339  # rad / s
    miu = 5.01
    A0_adiab = sqrt(miu) * Beta / pi  # [Hz]
    A0 = 1.5 * A0_adiab

    BW = Beta * miu
    N = int(round(duration / 1e-6))

    t = np.arange(1, N + 1) * system.rf_raster_time
    tt = t - (duration * center_pos)

    signal_tmp = complex(1, -1 * miu)
    signal_tmp2 = 1 / np.cosh(np.dot(Beta, tt))

    signal0 = np.power(signal_tmp2, signal_tmp)  # HyperSec    function
    signal_normal = signal0
    signal = A0 * signal_normal  # *1524

    rf_inv = SimpleNamespace()
    rf_inv.type = 'rf'
    rf_inv.signal = signal
    rf_inv.t = t
    rf_inv.freq_offset = freq_offset
    rf_inv.phase_offset = phase_offset
    rf_inv.dead_time = system.rf_dead_time
    rf_inv.ringdown_time = system.rf_ringdown_time
    rf_inv.delay = delay

    if rf_inv.dead_time > rf_inv.delay:
        rf_inv.delay = rf_inv.dead_time

    if rf_inv.ringdown_time > 0:
        t_fill = np.arange(1, round(rf_inv.ringdown_time / 1e-6) + 1) * 1e-6
        rf_inv.t = np.concatenate((rf_inv.t, rf_inv.t[-1] + t_fill))
        rf_inv.signal = np.concatenate((rf_inv.signal, np.zeros(len(t_fill))))

    negative_zero_indices = np.where(rf_inv.signal == -0.0)
    rf_inv.signal[negative_zero_indices] = 0

    return rf_inv


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # System limits
    system = Opts(max_grad=28, grad_unit='mT/m', max_slew=125,
                  slew_unit='T/m/s', rf_ringdown_time=20e-6,
                  rf_dead_time=100e-6, adc_dead_time=10e-6)

    seq = Sequence(system)

    # Sequence Parameters
    FOV = 200e-3  # [m]
    Ny = 128
    Nstartup = 11

    # Calculate TR
    _, TR, Ny_aq = bssfp_readout(Sequence(system), system, fov=FOV, Nstartup=Nstartup, Ny=Ny)

    # Calculate Crushers
    gzSpoil_INV = make_trapezoid('z', area=-3e3, duration=9.5e-3, system=system)
    FlatArea_crusher = [0, 0.4e3, 0.32e3, 0.27e3, 0.32e3, 0, 0.53e3, 0.32e3]
    FlatTime_crusher = [0, 5e-4, 6.3e-4, 6.3e-4, 6.3e-4, 0, 9.3e-4, 6.3e-4]
    crushers_in_between = 1

    # Create adiabatic inversion pulse
    rf_inv = make_hyperSec_pulse(system=system, duration=10.24 * 1e-3)
    InvDur = rf_inv.t[-1]  # duration

    Tdelay_trig = 550e-3  # [s]
    trig = make_trigger(delay=100e-6, duration=Tdelay_trig - (Nstartup + 1) * TR, system=system)
    trig_BetweenInversion = make_trigger(delay=100e-6, duration=Tdelay_trig, system=system)
    # Next pypulseq version add channel
    #trig = make_trigger(channel='physio1', delay=100e-6, duration=Tdelay_trig - (Nstartup + 1) * TR, system=system)
    #trig_BetweenInversion = make_trigger(channel='physio1', delay=100e-6, duration=Tdelay_trig, system=system)

    # Scheme 5(3)3
    NImage_perInv = [5, 3]
    Nrecover_btwInv = [3]
    TI_Vector = [0.100, 0.100 + 0.08]  # [s]     TI1 minimum TI of 100 ms, TI increment of 80 sec, Messroghli 2007
    TI_Vector_real = np.zeros(np.shape(TI_Vector))  # Initialize vector for final TI
    N_inversion = np.shape(NImage_perInv)  # Number of inversions

    print('Number of inversion:', N_inversion[0])
    print('Number of images per inversion:', NImage_perInv)
    print('Number of recovery heart beats between inversion:', Nrecover_btwInv)

    # Minimum Inversion Time
    TI_min_attainable = (Nstartup + Ny / 2) * TR + calc_duration(gzSpoil_INV)
    print('Minimum TI: %3.0f' % (TI_min_attainable * 1e3), 'ms')

    # Iterate between Inversions
    for nInv in np.arange(N_inversion[0]):
        print('----------- Inversion %3.0f -----------' % (nInv + 1))
        if nInv > 0:
            for nRec in np.arange(Nrecover_btwInv[nInv - 1]):
                seq.add_block(trig_BetweenInversion)

        # ***************************************************************************************************
        #                          Verify minimum possible TI
        #                   (nearest to TI_Vector[0], i.e. 100 ms)
        try:
            assert (all(TI_Vector[nInv] >= TI_min_attainable))  # verify if TI_vector value is possbile
            delayINV = TI_Vector[nInv] - TI_min_attainable  # calculate delay for TI
            TI_Vector_real[nInv] = TI_Vector[nInv]
            seq.addBlock(make_delay(delayINV))

        except:
            delayINV = 0
            if nInv == 0:
                TI_Vector_real[nInv] = TI_min_attainable
            else:
                TI_Vector_real[nInv] = TI_min_attainable + np.diff(TI_Vector)  # add increment i.e. 80 ms

            print('TI changed from ', TI_Vector[nInv] * 1e3, ' ms to %3.0f ms' % (TI_Vector_real[nInv] * 1e3))
        # ****************************************************************************************************

        trig_inv = make_trigger(delay=100e-6,
                                duration=Tdelay_trig - delayINV - InvDur - calc_duration(gzSpoil_INV) - (
                                        Nstartup + 1) * TR,
                                system=system)  # next pypulseq version add channel='physio1',
        seq.add_block(trig_inv)

        ##########################################################################################
        # Inversion pulse
        seq.add_block(rf_inv)
        ##########################################################################################

        seq.add_block(gzSpoil_INV)  # crusher

        if delayINV > 0:
            seq.add_block(make_delay(delayINV))  # Delay for min attainable TI

        # Add readouts for each inverion
        for nACQ in np.arange(NImage_perInv[nInv]):
            if nACQ > 0:
                gzCrusher_im = make_trapezoid('z', flat_area=FlatArea_crusher[nACQ], flat_time=FlatTime_crusher[nACQ],
                                              system=system)

                ##########################################################################################
                # Add trigger between images of the same inversion
                trig_crush = make_trigger(delay=100e-6,
                                          duration=Tdelay_trig - (Nstartup + 1) * TR - calc_duration(gzCrusher_im),
                                          system=system)  # next pypulseq version add channel='physio1',

                seq.add_block(trig_crush)
                ##########################################################################################

                if crushers_in_between:
                    seq.add_block(gzCrusher_im)  # spoiler before image acquisition
                else:
                    seq.add_block(make_delay(calc_duration(gzCrusher_im)))

            ##########################################################################################
            # bSSFP readout
            seq, TR, _ = bssfp_readout(seq, system, fov=FOV, Nstartup=Nstartup, Ny=Ny)
            ##########################################################################################

    print('---------------  Done  ---------------')
    print('First Inversion Time', TI_Vector_real, 's')

    #Plot complete sequence
    seq.plot()

    seq.write('pyProMyoT1.seq')