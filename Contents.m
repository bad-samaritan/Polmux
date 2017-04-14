% Optical simulation toolbox
% Version 0.1
%
% Author: Paolo Serena, University of Parma, 2009.
%
%
%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2009  Paolo Serena, <serena@tlc.unipr.it>
%			 
%    Optilux is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    Optilux is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%   reset_all         - Reset all global variables and initialize the simulation.
%   create_field      - Create the electric field.
%   pattern           - Create the sequence pattern with rules.
%   pat_encoder       - Symbols encoder.
%   pat_decoder       - Symbols decoder.
%   pat2stars	      - Convert an M-ary pattern into a complex constellation
%   samp2pat	      - Convert received samples into a pattern
%   stars2pat	      - Convert a complex constellation into a pattern
%   electricsource    - Create the electric modulating signal.
%   lasersource       - Multichannel laser transmitter.
%   linear_modulator  - Modulate the optical field with a linear modulator.
%   mz_modulator      - Modulate the optical field with a Mach-Zehnder Interferometer.
%   phase_modulator   - Modulate the optical field with a phase modulator.
%   qi_modulator      - Modulate the optical field using a QI Mach-Zehnder modulator
%   lpfilter          - Filtering with a lowpass filter.
%   fiber             - Optical fiber in the nonlinear regime.
%   fibergui          - Optical fiber in the nonlinear regime (GUI tool).
%   ampliflat         - Ideal Optical amplifier with ASE noise.
%   optfilter         - Optical filter.		  - 
%   receiver_ook      - Complete OOK receiver (POST fiber+OBPF+photodiode+LPF).
%   receiver_dpsk     - Complete DPSK receiver. (POST fiber+OBPF+MZ+LPF).
%   receiver_dqpsk    - Complete DQPSK receiver. (POST fiber+OBPF+MZs+LPF).
%   receiver_cohmix   - Complete COHerent MIXer receiver. (POST fiber+OBPF+MIX+PD+LPF). 
%   eval_eye          - Evaluate the eye opening for a non -coherent transmission.
%   eval_polar	      - Evaluate the polar diagram for a coherent transmission.
%   ber_kl            - Evaluate the BER for noncoherent transmission by Karhunen-Loeve method.
%   best_eye          - Search algorithm for the best eye opening.
%   best_sp           - Search algorithm for the best OSNR penalty vs. back-to-back.
%   ber_estimate      - Bit -error rate estimate by Monte Carlo simulation.
%   mc_estimate       - Monte Carlo estimation of a random variable mean and variance.
%   cmaadaptivefilter - Polarization demultiplexing filter using CMA algorithm
%   dsp4cohdec	      - Digital signal processing for a coherent receiver
%   easiadaptivefilter- Source separation filter using EASI algorithm
%   nmod              - N -modulus of an integer.
%   pow2phi 	      - Convert power into nonlinear phase.
%   phi2pow           - Convert nonlinear phase into power.
%   avg_power         - Evaluate the average power per symbol
%   corrdelay	      - System delay by cross -correlation measurement.
%   evaldelay         - Evaluate the group -delay of the filter. 
%   myfilter          - Filter device in the frequency domain.
%   pol_scrambler     - Rotates the SOP of signal samples on the Poincare sphere.
%   dop_meter	      - Computes the Degree Of Polarization of the Optical field.
%   inverse_pmd       - Inverse PMD matrix.
%   set_sop 	      - Sets the average State Of Polarization of a signal.
%   polarizer         - Linear optical polarizer
%   plotfield         - Plot the optical field.
%   plotfile          - Plot file from disk.
%   printfield        - Print the optical field to file.
%   ber2q	      - Convert bit -error rate in Q factor.
%   fprintfmsg        - Write a message into the file simul_out.
%   fastexp 	      - Calculates exp( i * x ) quickly.
%   saddle            - Evaluate the MGF saddle point.
%   comp_mex	      - compile all .c files into the directory.
%   fastshift         - Fast but simplified circular shift.
%   mdoc	      - Display Optilux HTML documentation in the browser
%   checkfields       - Check for valid input fields
