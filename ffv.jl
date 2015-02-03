# implementation of Fundamental Frequency Variation
# Based off the the code and papers of Kormel Lakowski
# Author : Boon Pang Lim
#   Date : Tue Feb  3 17:26:22 SGT 2015
# verified to be numerically similar to FFV.

# NB: this doesn't exactly match the code from Kormel -- fft computations seem off.

using WAV;

# Default settings and constants from ffv-1.0 code 
#
# tint           : 0.011
# text           : 0.009
# tsep           : 0.014
# fs             : 16000
# Tint           : 176
# Text           : 144
# Tsep           : 224
# winShape       : 2
# Nh             : 512
# NhPow2         : 512
# Nfft           : 1024
# Ng             : 512
# Nf             : 7
# tsepRef        : 0.008
# dCorrXFormType : 4
# nOutput        : 7

default_tint = 0.011;
default_text = 0.009;
default_tsep = 0.014;
default_tsepref = 0.008
default_fs   = 16000.0;
default_tfra = 0.008;
default_NhPow2 = 512;
default_Ng = 512;

# parameters for various window shapes and types
type WINPARAM
    bint :: Float64;
    mint :: Float64;
    bext :: Float64;
    mext :: Float64;
end

hamming_hamming = WINPARAM(0.54,0.46,0.54,0.46);
hanning_hanning = WINPARAM(0.5,0.5,0.5,0.5);
hamming_hanning = WINPARAM(0.5,0.5,0.54,0.46);

# creates a pair of windowing functions 
# we default to the hamming/hanning window. This splices a hamming half window
# on the extremity and a hanning window towards the centre.
function create_winpair(fs=16000,size=512;
    tint=default_tint,
    text=default_text,
    tsep=default_tsep,windows=hamming_hanning)

    # computes left and right windows for ffv analysis
    hL=zeros(size,1);
    hR=zeros(size,1);

    # scale by sampling rate
    Tint = tint * fs;
    Text = text * fs;
    Tsep = tsep * fs;
    TsepHalf=int(Tsep/2)

    # left half ranges
    t_r =   -(Text + TsepHalf):((-TsepHalf)-1);
    i_r = 1+t_r + (Text + TsepHalf);
    hL[i_r] = windows.bext + windows.mext* cos(pi*(t_r+TsepHalf)/ Text)

    t_r =   (-TsepHalf):((-TsepHalf+Tint)-1);
    i_r = 1+t_r + (Text + TsepHalf);
    hL[i_r] = windows.bint + windows.mint* cos(pi*(t_r+TsepHalf)/ Tint)

    # right half ranges
    t_r =  (TsepHalf-Tint):(TsepHalf-1);
    i_r = 1+t_r + (Text + TsepHalf);
    hR[i_r] = windows.bint + windows.mint* cos(pi*(t_r-TsepHalf)/ Tint)

    t_r =  (TsepHalf):(TsepHalf+Text-1);
    i_r = 1+t_r + (Text + TsepHalf);
    hR[i_r] = windows.bext + windows.mext* cos(pi*(t_r-TsepHalf)/ Text)

    return hL,hR;
end

# Linearly interpolates the magnitude spectrum for a fractional FFT bin index
# inputs:     sck - the fractional bin
#         fft_mag - magnitude spectrum
function interp_mag(sck,fft_mag)
    sck=abs(sck)
    sck1=floor(sck);
    sck2=ceil(sck);
    if sck1==sck2 then
        return fft_mag[1+sck1]
    else
        coeff1=abs(sck-sck2);
        coeff2=abs(sck-sck1);
        return coeff1*fft_mag[1+sck1]+coeff2*fft_mag[1+sck2]
    end
end

function interp_mag2(sck,fft_mag;nFFT=512)
    "linearly interpolate magnitued fft"
    if sck<0 then
        d_sck = floor(sck) + nFFT + 1
        p_sck = ceil(sck)  + nFFT + 1
        proxCoeff = abs(floor(sck) - sck)
    else
        d_sck = ceil(sck)  +1
        p_sck = floor(sck) +1
        proxCoeff = abs(ceil(sck) - sck)
    end
    distCoeff = 1.0-proxCoeff

    if p_sck > nFFT
        p_sck = 1
    end

    prox_mag = abs(fft_mag[p_sck])
    dist_mag = abs(fft_mag[d_sck])

    mag = proxCoeff*prox_mag + distCoeff*dist_mag
    return mag
end

# Computes FFV spectrum for an already windowed frame.
function ffvSpectralSlice(cwin,hL,hR,tsep=default_tsep,tsepRef=default_tsepref,Ng=default_Ng,NhPow2=default_NhPow2)

    # compute left and right windows
    lwin = hL.* cwin;
    rwin = hR.* cwin;

    # compute magnitude spectrum
    nFFT=length(hL);

    lfft_mag = abs(fft(lwin)[1:nFFT])
    rfft_mag = abs(fft(rwin)[1:nFFT])

    tsepRatio = tsep/tsepRef;

    s=zeros(Ng,1);
    # compute for descending pitch - stretch right frame spectrum to match left
     
    #for r in -Ng/2:-1
    for r in -Ng/2:-1
        rho=2.0^ (4.0*(-(abs(r)/Ng)*tsepRatio));
        product=0; normL=0; normR=0
        for k in -NhPow2/2:(NhPow2/2-1)
            sck=rho*k
            fL_mag=abs(lfft_mag[abs(k)+1])
            fR_mag=interp_mag(sck,rfft_mag)
            product += fL_mag*fR_mag;
            normL += fL_mag*fL_mag;
            normR += fR_mag*fR_mag;
        end
        s[r+Ng/2+1] = product/sqrt(normL*normR+1e-9);
    end

    # compute for ascending pitch - stretch left frame spectrum to match right
    for r in 0:Ng/2-1
        rho=2.0^ (4.0*(-(abs(r)/Ng)*tsepRatio));
        product = 0; normL = 0; normR = 0;
        for k in -(NhPow2/2):(NhPow2/2-1)
            sck = rho * k;      
            fL_mag = interp_mag(sck,lfft_mag)
            fR_mag = rfft_mag[abs(k)+1];
            product += fL_mag*fR_mag;
            normL += fL_mag*fL_mag;
            normR += fR_mag*fR_mag;
        end
        s[r+Ng/2+1] = product/sqrt(normL*normR+1e-9);
    end
    return s;
end

# Computes FFV spectrum for an audio signal
# options - tfra - frame step in seconds
function ffvSpectrum(fs,samples;tfra=default_tfra)
    hL,hR=create_winpair();

    fr_step = int((tfra*fs));
    fr_size=512;
    nfr = length(1:fr_step:length(samples))
    
    # perform preemphasis and zero pad the signal 
    nsamples=zeros((nfr+10)*fr_step,1);
    nsamples[2:length(samples)] = -0.97*samples[1:length(samples)-1]+samples[2:length(samples)]
    nsamples=nsamples*32768; # correct for data conversion from short to float

    out_spectrum=zeros(512,nfr);

    Ng = 512;   # defaults for number of vanishing points

    for k in 1:nfr
        i=1+(k-1)*fr_step
        cwin=nsamples[i:(i+fr_size-1)]
        out_spectrum[:,k] = ffvSpectralSlice(cwin,hL,hR);
    end
    return out_spectrum;
end

# currently only works for default 7 point filter bank and 512 point FFV
function ffvCreateFreqFilterbank(nf=7,nFFV=512)
    # for 512-bin FFV spectra
    # 117 to 139
    # 245 - 251
    # 249 - 255 (?!? Bug?)
    # 254 - 258
    # 257 - 263
    # 261 - 267
    # 373 to 396

    F={ [ 118:140 ], [ 246:252 ], [ 250:256 ],[ 255:259 ], [ 258:264 ],[ 262:268 ], [ 374:397 ] };

    function rect(arr)
        y=ones(length(arr),1)/float(length(arr))
    end

    function trapezoid(arr)
        y=ones(length(arr),1);
        y[1] = 0.5;
        y[length(y)] = 0.5;
        y=y./float(length(arr)-1);
    end

    G=cell(size(F))
    G[1]   = [F[1] rect(F[1])]
    G[2:6] = map(x->[x trapezoid(x)],F[2:6])
    G[7]   = [F[7] rect(F[7])]
    return G;
end

# Computes FFV for an audio signal, applying 7-point filterbank to the spectrum
# options - tfra - frame step in seconds
function ffv(fs,samples;tfra=default_tfra)
    Fb = ffvCreateFreqFilterbank()
    S=ffvSpectrum(fs,samples;tfra=tfra);
    print (size(S))
    out=zeros(7,size(S,2))
    for i in 1:length(Fb)
        X=Fb[i]
        out[i,:] = sum(S[X[:,1],:].*X[:,2],1)
    end
    return out;
end
