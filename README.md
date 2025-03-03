# adaptive_unbias_lateral_filter
oct image denoise
Enhancement and Bias Removal of Multiframe Optical

Coherence Tomography Images: an Iterative Approach

via Adaptive Bilateral Filtering

PV Sudeepa,∗

, S Issac Niwasb,c, P Palanisamya

, Jeny Rajand

,

Y. Xiaojunb,c, X. Wangb,c, Y. Luob,c, L. Liuc,e

aDepartment of Electronics and Communication Engineering, National Institute of

Technology - Tiruchirappalli, Tamil Nadu, India

bSchool of Electrical and Electronic Engineering, Nanyang Technological University

(NTU), Singapore 639798, Singapore.

cCentre for Optical Fibre Technology (COFT), The Photonics Institute (TPI), Nanyang

Technological University (NTU), Singapore 639798, Singapore.

dDepartment of Computer Science and Engineering, National Institute of Technology

Karnataka, Surathkal, India

eSchool of Chemical and Biomedical Engineering, Nanyang Technological University

(NTU), Singapore 639798, Singapore.

Abstract

Optical coherence tomography (OCT) has continually evolved and expanded

as one of the most valuable routine tests in ophthalmology. However, noise

(speckle) in the acquired images causes visual quality degradation of OCT

images and makes it difficult to distinguish structural details in the acquired

image. In this article, an iterative approach based on bilateral filtering is

proposed for speckle reduction in multiframe OCT data. Gamma noise model

is assumed for the observed OCT image. First, the adaptive version of the

conventional bilateral filter is applied to ameliotate the multiframe OCT data

and then the bias due to noise is reduced from each of the filtered frames.

∗Corresponding author

URL: spvnitt@gmail.com (PV Sudeep)

Preprint submitted to Computers in Biology and Medicine February 1, 2016

These unbiased, filtered single frames are then refined using an iterative

approach. Finally, these refined single frames are averaged to produce the

enhanced OCT image. Experimental results on phantom images and real

OCT retinal images demonstrate the effectiveness of the proposed filter.

Keywords: Bilateral filtering, denoising, maximum likelihood estimation,

optical coherence tomography, speckle noise.

1. Introduction

Optical coherence tomography (OCT) is a noninvasive imaging modality

that yields high resolution images of tissue structures and cross-sectional

imaging of many biological systems [1]. As an optical equivalent of ultrasound

(US) imaging, it uses echoes of light rather than sound [2]. Major OCT

applications include detection of eye diseases and skin disorders [3, 4].

In ophthalmology, assessment and observation of visual ailments require

details about the inner retinal structure of the eye. OCT images are used

to visualize the morphological structures of the retina and iris macula and

evaluate the retinal nerve fiber layer thickness in the peripapillary region

as well as other inner tissue lines associated with the anterior and poste￾rior sections of the eye [5]. The ophthalmologists can observe and measure

the anatomical correspondence of the intra-retinal layers viz. retinal nerve

fibre layer (RNFL), ganglion cell layer (GCL), inner-plexiform layer (IPL),

inner-nuclear layer (INL), outer-plexiform layer (OPL), inner and outer seg￾ment of the photoreceptors (IS/OS), the retinal pigment epithelium (RPE)

and choriocapillaris (CC)) by performing anterior and posterior imaging of

the eye using OCT [6, 7]. These are clinically important in case of many

2

ophthalmic diseases such as early detection and diagnosis of eye diseases,

including glaucoma, age-related macular degeneration, macular edema, and

diabetic retinopathy.

Nevertheless, the visual quality of OCT images can be severely damaged

by speckle noise arising due to the interference of photons that undergo

multiple scattering in reverse and forward direction, when photons propagate

inside tissue [8, 9]. It is crucial to remove unwanted speckle noise from

OCT images since the automated segmentation of the inner tissue lining

of eye images (e.g., optic nerve head layers, retinal layers, drusens) and its

measurement are clinically important for the diagnosis of eye diseases [10].

Apart from OCT images, speckle noise reduction is a significant research

topic in other imaging modalities such as synthetic aperture radar (SAR)

and US images.

In literature, various speckle reduction algorithms have been explored for

OCT images. In [11], Rogowska et. al have discussed techniques such as

image averaging, mean, median and Gaussian filters for removing unwanted

components from OCT images. The Lee filter [12], Frost filter [13] and Kuan

filter [14] are the most widely discussed spatial adaptive filters to attenuate

speckle noise. Enhanced versions of the Lee and Frost filters are proposed in

[15]. Computationally efficient speckle noise reduction algorithms using an

adaptive Wiener filter [16] and diffusion filtering techniques[17, 18, 19, 20]

have also been proposed speckle reduction.

In the past several methods are proposed for OCT denoising. A complex

diffusion filter was proposed in [21] and its improved adaptive technique was

discussed in [22]. A combined filter by integrating the benefits of PDE-based

3

approach and the wavelet transform is employed in [23]. Multi-resolution

based spatially adaptive discrete wavelet filter [24], complex wavelet filter

[25], and curvelets filters [26] are other significant methods proposed in the

literature to remove speckle noise. A. Ozcan et.al have compared various

filters including non-orthogonal wavelet transform based filters together with

enhanced Lee and adaptive Wiener filters for OCT speckle denoising [27]. I￾divergence regularization approach [28] and general Bayesian estimations [29]

have also been proposed for speckle reduction in OCT imaging. Moreover,

different patch based non-local recovery paradigms have recently introduced

for speckle noise reduction in OCT images [30, 31, 32].

A better statistical characterization of the observed (noisy) images can

help in developing effective despeckling methods. In the literature, different

statistical models such as Gaussian, Rayleigh and Log-normal have been in￾vestigated to address the restoration of speckled images [11, 33]. In this work,

a Gamma model is followed for multi-frame OCT images and a three parame￾ter Gamma distribution function is used to fit the observed OCT data. Three

parameter Gamma distribution is able to include several well known models

such as exponential, Rayleigh, Gamma, Chi-square and normal distribution

as subfamilies.

In this paper, we propose an adaptive and unbiased bilateral (AUB) filter

to account for the Gamma characteristics of the data and then, incorporated

the idea of iterative filtering to refine the AUB filter result. This iterative

AUB (IAUB) filtering technique is performed on each single frame, and the

outputs of the IAUB filter are averaged to obtain the final despeckled image.

The excellent functioning of the proposed filter is well validated by experi-

4

ments using both simulated and real OCT images and the simulations show

that it provides better speckle removal without affecting noticeable structures

in the image.

This paper is organized as follows: The noise model for OCT images is

briefly described in Section 2, and the proposed methodology is presented in

Section 3. Section 4 discusses the experimental results and the evaluation

process. Finally, we conclude the paper in Section 5.

2. Noise Characteristics in OCT images

The quality of images obtained by any image acquisition system employ￾ing coherent wave fields is degraded by speckle noise, an insidious form of

noise produced by the interference of waves with random phases [34, 35]. The

speckle noise manifests itself as a fast fluctuation of the detected intensity (or

field envelope) over the spatial extent of the image, conveying a granular tex￾ture [36]. High-speed OCT imaging quality is restricted due to the presence

of speckle noise [37]. With a large number of polarized quasimonochromatic

waves with random phase, a fully developed speckle pattern is formed [38].

In [39], Goodman et. al. find a relation for the speckle field using a

Rayleigh distribution function. More attempts to model speckle patterns are

found in the US literature than in the OCT literature. Despite the fact that

a multiplicative model and empirical distributions with useful probabilistic

models, e.g., Nakagami, generalized Nakagami, Weibull and Rician inverse

Gaussian (RiIG) distributions, have been studied [40, 41, 42, 43, 44, 45] for

modelling speckle noise in ultrasound, Vegas-S´anchez-Ferrero et. al. [46]

shows that Gamma distribution, which is a good approximation for the

5

weighted sum of Rayleigh variables, accurately fits interpolated output B￾mode US images, when compression-less data with fully developed speckle is

considered.

Motivated by the above studies, we assume the noise field in OCT images

is due to fully developed speckle (as discussed in [34]), and hence, we model

the OCT data distribution by Gamma probability density function (PDF).

Let us consider that there are k frames in the OCT data, and the dimension

of each frame is r ×c. If the noisy OCT data, their noiseless version, and the

noise are denoted as M, C and N, respectively, we can express the relation

between them as in Eq. (1),

M = C + N (1)

where the dimensions of the vectors M, C, and N are all (r × c) × k. In

other words, N denotes a fading variable and takes random numbers from the

Gamma distribution with shape and scale parameters ρ and β respectively.

3. Methodology

In this section, we concisely present the concept of bilateral filtering,

and then propose two filters to reduce speckle from the Gamma distributed

OCT data, namely, an adaptive and unbiased bilateral (AUB) filter and the

iterative adaptive and unbiased bilateral (IAUB) filter. These filters are

suitable for enhancing both single-frame and multiframe OCT images.

The conventional bilateral (CB) filter is a nonlinear spatial filter devel￾oped by Tomasi and Manduchi [47]. CB filter is very popular for edge pre￾served smoothing of images, which are Gaussian distributed. The kernel of

6

the CB filter is composed of two components, namely a range filter kernel

and domain filter kernel. The response of the CB filter at a pixel location m

can be calculated as [47]

I

ˆ(m) =

Z

1 X

n∈N(m)

wD (m, n)wR (m, n) I(n), (2)

where N(m) represents the neighbourhood region around m , and n is the

position in the neighbourhood. The normalization constant Z is given by

Eq. (3),

Z =

X

n∈N(m)

wD (m, n)wR (m, n). (3)

The weight function wD is linked to the domain filter and measures the

photometric similarity of a pixel around its neighbourhood. Similarly, the

other weight function wR is related to range filtering and the computed

weights are proportional to the radiometric distance around the neighbour￾hood of a pixel.These weight functions can be defined as in Eq. (4) and in

Eq. (5),

wD (m, n) = exp −

|m

2

−

σd

2

n|

2

!

(4)

and

wR (m, n) = exp −

|Im

2

−

σr

2

In|

2

!

(5)

where Im and In are the intensities at m and n, respectively. The geometrical

spread σd in the domain determines the amount of low pass filtering in the

form of blur, and its optimal value is relatively insensitive to noise. However,

the photometric spread σr in the image range is almost linearly related to

the true noise standard deviation σ [48, 49] .

7

3.1. Adaptive and unbiased bilateral (AUB) filter

The first step to derive the adaptive and unbiased bilateral (AUB) filter is

the conversion of the CB filter to an adaptive bilateral filter (ABF) by choos￾ing the filter parameters in an adaptive fashion. Even though two spreads

determine the performance of the CB filter, we need to make the photomet￾ric spread adaptive to reduce speckle because the other spread (geometric)

is insensitive to noise. Hence, we define σr as in Eq. (6):

σr = τ × σˆ (6)

where [50]

ˆσ

2 = mode n σ

ˆ2

i,jo (7)

σ

ˆ2 denotes the estimated noise variance from the given input image. The

preferred value for σd is in between 1.2 and 2.1 [48]. In our experiments to

find the value of τ , we fixed the σd value to 1.7, and experimentally find that

the optimum value for τ is 0.01. Now, we modify the ABF filter by reducing

the bias introduced due to noise in the image and propose the AUB filter.

In order to get an accurate estimate of the bias, we have developed an ML

estimation based method relying on the Gamma distributed nature of the

observed OCT data. In section 3.2, we will discuss the ML estimation of the

bias term.

3.2. Bias calculation using the maximum likelihood estimation method

We assume that the noisy OCT image is Gamma distributed, as in Eq. (1),

and hence, we model it with the three parameter Gamma distribution de￾scribed by Cohen and Whitten in [51]. Fig. 1 displays the distribution of a

8

homogeneous region of a real swine OCT image (shown in Fig. 5 (a)). Com￾parison with the true Gamma distribution with parameters estimated from

the selected homogeneous region shows that the data distribution can be

modeled with a Gamma PDF. Let u = u1, u2, ..., uk be k statistically inde￾pendent random observations from a Gamma distribution. Then, the PDF

of the observations is given by[51]

P(u|γ, ρ, β) =







(u−γ)

ρ−1

Γ(ρ)βρ exp  −

(u−

β

γ)



∀γ < u < ∞,

ρ > 0, β > 0

0 otherwise

(8)

where Γ (·) denotes the Gamma function and γ, ρ and β represent the loca￾tion, shape and scale parameters respectively.

The mean µu and variance σu of the above mentioned distribution are

µu = γ + ρβ (9)

σu

2 = ρβ2

(10)

respectively. Note that, the mean µu of the observed signal u has two com￾ponents, in which γ corresponds to the mean value of the noiseless signal and

the product term ρβ represents the mean of the Gamma distributed pixels

in the background, i.e., pixel locations in which the signal intensity is zero.

In other words, the product term ρβ represents a shift in the mean value of

the noiseless signal, and in this article, we have called this term the bias.

In order to compute the bias, the estimation of the distribution param￾eters, ρ and β is essential and can be estimated by maximizing the corre￾sponding likelihood function L , or equivalently log L with respect to ρ and

9

β. Since the parameter γ represents the noiseless signal component, it is un￾known at the stage of the bias computation. Hence, Eq. (9) can be rewritten

as,

γ = µu − ρβ (11)

Then, the unknown noiseless image parameter γ in Eq. (8) can be replaced

by Eq. (11).

Under the condition of γ < u < ∞, ρ > 0, β > 0, the ML estimator can

be defined as follows:

{ρML, βML} = arg max

ρ,β

(log L) (12)

where

log L =

nX

i=1

log P(ui

|ρ, β)

= (ρ − 1)

nX

i=1

log (ui − µu + ρβ)

−

β

1

nX

i=1

(ui − µu + ρβ) (13)

−n [ log Γ (ρ) + ρ log β ]

and ρML and βML are the estimated shape and scale parameters, respectively.

These parameters have fixed values throughout all regions in the input

image, and hence, we can compute them from any region in which the un￾derlying intensity is constant. If we assume that there exist many such

piecewise homogeneous regions then ρ and β can be estimated as ˆρML =

mode ￾ ρML(i,j)

	 and βˆML = mode ￾ βML(i,j)

	

, where ˆρML and βˆML are

the estimated values of ρ and β for the given image. Fig. 2 shows the distri￾bution of local estimates of ρ and β from a real swine OCT image (shown in

10

Fig. 5 (a)). Now, the bias(B) can be written as

B = ˆρMLβˆML (14)

This bias increases as the amount of noise in the OCT image increases. The

bias computation for each image is tedious and time consuming in a multi￾frame OCT dataset because the process has to be repeated for all the frames.

Since the experimental set up is same, bias will be more or less the same in

all frames. In any case, averaging bias from randomly selected frames will

help to reduce the estimation errors associated with the ML method.

It can be observed from Eq. (9) and Eq. (14) that the mean of the observed

Gamma distributed OCT image is biased by the speckle noise and hence, bias

must be subtracted from the mean filtered image to achieve the unbiased,

filtered image.

For each pixel of the discrete noisy image I =

 I(n)|n ∈ R

N

	 , the pro￾posed AUB estimator can be defined as follows:

I

ˆ

AUB (I (n)) = max  I

ˆ

ABF (I (n)) − B, 0

 (15)

Therefore, for a k-frame (multiframe) OCT data, the final averaged AUB

filter output is given by,

I

ˆ

avg =

k

1 X

k

I

ˆ

AUB,k (Ik (n))

=

k

1 X

k

max  I

ˆ

ABF,k (Ik (n)) − B, 0

 (16)

11

Algorithm 1 The IAUB filtering method for speckle reduction

1: Input: noisy multi-frame OCT images (say, number of frames as R)

2: Compute bias B by using Eq. (14).

3: for p=1 to R do

4: Estimate ˆσ from p

th OCT image I

p by using Eq. (7).

5: Calculate σr

p with Eq. (6).

6: Smooth I

p by the CB filter defined in Eq. (2) to obtain I

ˆp

.

7: Find the AUB filter output I

ˆ

1

p

from I

ˆp by using Eq. (15).

8: Set I

p = I

ˆ

1

p

.

9: for q=2 to iter do

10: Estimate the new value for ˆσ and hence, σr

p by using

Eq. (6) and Eq. (7) from I

p

.

11: Find the refined image I

ˆ

q

p

using the CB filter defined in Eq. (2).

12: Set I

p = I

ˆ

q

p

13: end for

14: IIAUB,p = I

p

.

15: end for

16: Compute IIAUB = ˆIavg = p

1 P

p

IIAUB,p.

17: Output: IIAUB ⇐ denoised output image.

3.3. Iterative adaptive and unbiased bilateral (IAUB) filtering of multi-frame

OCT images

The quality of the image denoised with the proposed AUB can be fur￾ther improved by applying the bilateral filter iteratively. Fig. 3 displays the

distribution of a homogeneous region of the AUB filter output OCT image

shown in Fig. 5 (f)). Comparison with the true Gaussian distribution with

12

parameters estimated from the selected homogeneous region shows that the

data distribution can be modeled with a Guassian PDF. As a well known

filter to reduce Gaussian noise, the bilateral filter can be used to refine the

result further. So, an iterative approach is introduced in this paper to refine

the AUB filter output and the noise variance is adaptively adjusted.The pro￾posed iterative adaptive and unbiased bilateral (IAUB) filtering methodology

is listed in Algorithm 1.

4. Experimental Results and Discussions

In this section, we report the experimental results obtained with synthetic

phantom and real OCT images. The results were evaluated through visual

inspection and quantitative analysis, and the proposed filter is compared

with other state-of-the-art methods.

The various parameters used for the filters discussed in this article are as

follows:

1. Bayesian estimation [29]: σspatial=1.5; window-size parameter =3; num￾ber of samples used = 70; Sigma factor = 2;

2. Wavelet multi-frame filter [37]: Window size for calculating correlation

= 5 × 5; wavelet basis = Haar; maximum decomposition level = 1;

weightmode-parameter = 1, [K P R] = [1 5 2].

3. Wiener filter [16]: Neighborhoods size (to estimate the local image

mean and standard deviation) = [5, 5].

4. Bilateral filter [47]: Half width of the filter=5; geometrical spread

σd=1.7; photometric spread σr=0.1.

13

5. Proposed AUB filter: Half width of the filter=5; geometrical spread

σd=1.7; photometric spread σr = 0.01 × σ.

6. Proposed IAUB filter: Half width of the filter=5; geometrical spread

σd=1.7; photometric spread σr = 0.01 × σ.

We found that the number of iterations can vary the performance of the

IAUB filter and the optimum values for number of iterations used, in our

experiments, with different bias B are given below

iter =







2 when

3 20 < B <

B < 20

40

,

,

4 B > 40.

(17)

4.1. Experiments on simulated images

This section reports the results of the qualitative and quantitative analysis

on a synthetic phantom image of size 256 × 256. Fig. 4(a) shows the Ground

Truth (GT) and Fig. 4 (b) shows the GT corrupted by Gamma noise with

parameters (ρ, β) = (4, 4). From Fig. 4, it can be observed that the image

denoised with the proposed IAUB filter (see Fig. 4 (n))is much closer to

the GT than the image filtered with other speckle removal techniques. For

better visual comparison, we have also shown corresponding residual images

in Fig. 4.

Qualitatively, the proposed IAUB filter provides better preservation of

fine structures, good contrast between different regions and better speckle

removal than the other methods under consideration. Our quantitative eval￾uation relies on objective measures such as the peak signal to noise ratio

(PSNR) [52], mean structural similarity index matrix (MSSIM) [53], Bhat￾tacharrya coefficient (BC) [54] and Pratts Figure of Merit (Pratt-FOM) [55].

14

For the quantitative analysis, the phantom image in Fig. 4 (a) is artificially

corrupted with Gamma random noise generated for different shape and scale

parameters, and hence, with various values of the bias B = 6 to 56. The

quantitative analysis based on PSNR, MSSIM, BC and FOM is summarized

in Table 1, Table 2, Table 3 and Table 4, respectively. It can be observed from

the tables that the proposed IAUB filter has a better performance over other

state-of-the-art methods, like Bayesian estimation, the wavelet multi-frame

filter, the Wiener filter, the wavelet filter, and the bilateral filter, along with

the proposed AUB filter.

4.2. Experiments on real OCT images

In order to validate the proposed methodology and to do a visual check

on the behaviour of the filter, experiments were carried out on real datasets

obtained by spectral domain OCT (SD-OCT) systems. In this section, we

present the results of two sample datasets: first, a real OCT data set of

a swine retina, which is acquired by our self-customized OCT system from

the Centre for Optical Fibre Technology (COFT), Nanyang Technological

University (NTU), Singapore, and the second, a pig’s eye dataset, which is

available online [37].

The first set of swine retinal images was obtained using an SD-OCT sys￾tem, which adopts a superluminescent diode (SLD) (Superlum Broadlighters

T-850-HP) with a center wavelength of λc = 850nm and a full width at

half maximum (FWHM) bandwidth of ∆λ = 165nm as the laser source. It

achieves a resolution of 2.3µm in air in the axial direction, and a resolution of

2.5µm in the transverse direction. The sensitivity of the system is measured

to be 105.6dB with 5.65mW light incident on the sample. The imaging speed

15

of the system is set to be 20 frames per second, with each frame consisting of

512 A-lines. We use a total of 20 images in our experiment for methodology

validation, and each image consists of 1024 transverse × 512 axial pixels,

covering an area of 0.495mm × 2.24mm.

The set of pig retinal images were obtained from the existing work [37].

This dataset was originally obtained from the Spectralis HRA & OCT (Hei￾delberg Engineering) device. The imaging speed of the system is 13 frames

per second, with each frame consisting of 768 A-scans. It achieves a reso￾lution of 7µm in air in the axial direction, and a resolution of 3.87µm in

the transverse direction. The dataset comprises of 35 sets, with 13 frames

included in each set sharing the same imaging position and hence, a total of

455 images in the dataset. Each image consists of 768 transverse × 496 axial

pixels with a resolution of 7µm in the axial direction. As stated in [37], the

reconstruction quality gradually improves using increasing number of frames

from 2 to 8. This trend is much less obvious with more than 8 frames. We

used total 13 frames (images) of a single set of data in the following experi￾ments to compare our approach with other methods.

The estimate of the shape and scale parameters for the pig eye dataset

and the swine eye retinal dataset are  ρ, ˆ βˆ

 = (2, 18) and  ρ, ˆ βˆ

 = (5, 8),

respectively. Hence, their corresponding bias B are 36 and 40, respectively.

The denoising results for the real datasets can be observed in Fig. 5 and Fig. 6.

Notice that, the images shown in Fig. 5 and Fig. 6 are obtained by averag￾ing the corresponding denoised frames. The despeckled image in Fig. 5(g)

and Fig. 6(g) shows the potential of our proposed IAUB filter in diminishing

the speckle noise associated with the ophthalmic OCT images. Its result

16

excels over the despeckled images by Bayesian estimation, the wavelet multi￾frame filter, the Wiener filter, the bilateral filter and our proposed AUB filter

(see Fig. 5(b) to Fig. 5(f), Fig. 6(b) to Fig. 6(f)). The retinal layer structures

(see Fig. 5(g) and Fig. 6(g)) are more visible using the proposed IAUB filter￾ing method and it has a very good edge preservation capability. It can be

observed that the results of the quantitative evaluation, shown in Table 1,

Table 2, Table 3 and Table 4, match well with the qualitative results in Fig. 5

and Fig. 6.

As an additional measure for evaluation, the expert ophthalmologists

scoring on the quality of the various filter outputs in the real OCT images

in Fig. 5 and Fig. 6 are provided in Table 5. Note that, the ophthalmologists

have graded the images on a scale of 1 − 10, for which 1 and 10 are used

to represent the poorest and best results, respectively. It can also be noted

through visual inspection that the image contrast is much better for our pro￾posed method when compared to the images denoised with other methods.

We examine how various filter parameters affect the performance of our

proposed method. The proposed filter has mainly four input parameters:

the shape and scale parameters ρ and β of the Gamma distribution are

computed in an automated way and the other parameters are the search

window size and similarity window size. Most importantly, the parameters ρ

and β contribute the bias and their estimation must be feasible. Otherwise,

the estimation errors in these parameters will considerably reduce the quality

of the filtered result. We applied the proposed IAUB technique to other OCT

images and obtained similarly good results.

One drawback of the proposed filters over other filters is its computational

17

burden. The AUB filter has a time complexity as that of the conventional

bilateral filter. However, additional time is required for computing the noise

level in the image (i.e., for the ML estimation of scale and shape parameters

of Gamma distribution).

5. Conclusion

In this paper, new speckle filters based on Gamma statistics have been

presented. First, we assume a Gamma distribution for the data in the noisy

frames of the given OCT data, and a local ML estimator is used to esti￾mate the Gamma distribution parameters, and hence the bias term due to

noise. In the proposed AUB filter, the calculated bias value is subtracted

from an adaptive version of the CB filter. The proposed IAUB filter has

an inherent AUB filter, and an iterative approach further refines the AUB

filter output. For this purpose, Gaussian distribution is assumed (based on

the experimental analysis) for the AUB filter output, and adaptive bilateral

filter is applied iteratively considering the superior performance of bilateral

filter on Gaussian distributed data. Experiments were carried out on both

simulated and real OCT images to validate the performance of the proposed

method. Evaluations based on the PSNR, mean SSIM, BC and Pratt-FOM

on synthetic images with various noise levels and visual inspections shows

that the proposed frameworks provide better performance than the state-of￾the-art methods. Also, results indicate that the proposed method achieves

considerable speckle noise reduction as well as better detail preservation for

the retinal images.

18

Acknowledgement

We sincerely appreciate funding support from Nanyang Technological

University, National Research Foundation Singapore (NRF2013NRF-POC001-

021 and NRF-CRP13-2014-05), National Medical Research Council Singa￾pore (NMRC/ CBRG/ 0036/ 2013), and Ministry of Education Singapore

(MOE2013-T2-2-107) for Nanyang Technological University affiliated au￾thors.

We would like to thank Dr. Victor Koh, MMed (Opthalmology), De￾partment of Ophthalmology, National University Health System (NUHS),

Singapore and Dr. S. Sujatha, MD, Professor of Ophthalmology, Joseph

Eye Hospital, Tiruchirappalli, India for providing their evaluation of various

methods on real OCT images used in this work.

References

[1] W. Drexler, J. G. Fujimoto, Optical Coherence Tomography: Technol￾ogy and Applications, Optical Coherence Tomography, Springer Inter￾national Publishing, 2015.

[2] J. Fujimoto, C. Pitris, S. Boppart, M. Brezinski, Optical Coherence To￾mography: An emerging technology for biomedical imaging and optical

biopsy, Neoplasia (New York, NY). 2 (1-2) (2000) 9–25.

[3] A. Fercher, C. Hitzenberger, G. Kamp, S. El Zaiat, Measurement of

intraocular distances by backscattering spectral interferometry, Optics

Communications 117 (12) (1995) 43 – 48.

19

[4] Y. Hori, Y. Yasuno, S. Sakai, M. Matsumoto, T. Sugawara, V. D. Mad￾jarova, M. Yamanari, S. Makita, T. Yasui, T. Araki, M. Itoh, T. Yata￾gai, Automatic characterization and segmentation of human skin using

three-dimensional optical coherence tomography, Opt. Express 14 (5)

(2006) 1862–1877.

[5] M. E. van Velthoven, D. J. Faber, F. D. Verbraak, T. G. van Leeuwen,

M. D. de Smet, Recent developments in optical coherence tomography

for imaging the retina, Progress in Retinal and Eye Research 26 (1)

(2007) 57 – 77.

[6] J. Mo, X. Yu, L. Liu, High Resolution Optical Coherence Tomography

for Bio-Imaging, in: M. Olivo, U. Dinish (Eds.), Frontiers in Biophoton￾ics for Translational Medicine, Vol. 3 of Progress in Optical Science and

Photonics, Springer Singapore, 2016, pp. 161–208.

[7] V. Guedes, J. S. Schuman, E. Hertzmark, G. Wollstein, A. Correnti,

R. Mancini, D. Lederer, S. Voskanian, L. Velazquez, H. M. Pakter,

T. Pedut-Kloizman, J. G. Fujimoto, C. Mattox, Optical coherence to￾mography measurement of macular and nerve fiber layer thickness in

normal and glaucomatous human eyes, Ophthalmology 110 (1) (2003)

177 – 189.

[8] F. K. Horn, C. Y. Mardin, R. Laemmer, D. Baleanu, A. M. Juenemann,

F. E. Kruse, R. P. Tornow, Correlation between Local Glaucomatous

Visual Field Defects and Loss of Nerve Fiber Layer Thickness Measured

with Polarimetry and Spectral Domain OCT, Investigative Ophthalmol￾ogy and Visual Science 50 (5) (2009) 1971.

20

[9] M. Szkulmowski, I. Gorczynska, D. Szlag, M. Sylwestrzak, A. Kowal￾czyk, M. Wojtkowski, Efficient reduction of speckle noise in optical co￾herence tomography, Opt. Express 20 (2) (2012) 1337–1359.

[10] P. P. Srinivasan, S. J. Heflin, J. A. Izatt, V. Y. Arshavsky, S. Farsiu, Au￾tomatic segmentation of up to ten layer boundaries in SD-OCT images

of the mouse retina with and without missing layers due to pathology,

Biomed. Opt. Express 5 (2) (2014) 348–365.

[11] J. Rogowska, M. E. Brezinski, Image processing techniques for noise re￾moval, enhancement and segmentation of cartilage OCT images, Physics

in Medicine and Biology 47 (4) (2002) 641.

[12] J.-S. Lee, Digital image enhancement and noise filtering by use of local

statistics, Pattern Analysis and Machine Intelligence, IEEE Transac￾tions on PAMI-2 (2) (1980) 165–168.

[13] V. Frost, J. Stiles, K. Shanmugan, J. Holtzman, A Model for Radar Im￾ages and Its Application to Adaptive Digital Filtering of Multiplicative

Noise, Pattern Analysis and Machine Intelligence, IEEE Transactions

on PAMI-4 (2) (1982) 157–166.

[14] D. Kuan, A. Sawchuk, T. Strand, P. Chavel, Adaptive noise smooth￾ing filter for images with signal-dependent noise, Pattern Analysis and

Machine Intelligence, IEEE Transactions on PAMI-7 (2) (1985) 165–177.

[15] T. Loupas, W. McDicken, P. Allan, An adaptive weighted median fil￾ter for speckle suppression in medical ultrasonic images, Circuits and

Systems, IEEE Transactions on 36 (1) (1989) 129–135.

21

[16] J. S. Lim, Two-dimensional Signal and Image Processing, Prentice-Hall,

Inc., Upper Saddle River, NJ, USA, 1990.

[17] P. Perona, J. Malik, Scale-space and edge detection using anisotropic

diffusion, Pattern Analysis and Machine Intelligence, IEEE Transactions

on 12 (7) (1990) 629–639.

[18] Y. Yu, S. Acton, Speckle reducing anisotropic diffusion, Image Process￾ing, IEEE Transactions on 11 (11) (2002) 1260–1270.

[19] S. Aja-Fernandez, C. Alberola-Lopez, On the estimation of the coeffi-

cient of variation for anisotropic diffusion speckle filtering, Image Pro￾cessing, IEEE Transactions on 15 (9) (2006) 2694–2701.

[20] K. Krissian, C. Westin, R. Kikinis, K. Vosburgh, Oriented Speckle Re￾ducing Anisotropic Diffusion, Image Processing, IEEE Transactions on

16 (5) (2007) 1412–1424.

[21] H. Salinas, D. Fernandez, Comparison of PDE-based nonlinear diffusion

approaches for image enhancement and denoising in Optical Coherence

Tomography, Medical Imaging, IEEE Transactions on 26 (6) (2007) 761–

771.

[22] R. Bernardes, C. Maduro, P. Serranho, A. Ara´ujo, S. Barbeiro,

J. Cunha-Vaz, Improved adaptive complex diffusion despeckling filter,

Opt. Express 18 (23) (2010) 24048–24059.

[23] A. Ogier, P. Hellier, C. Barillot, Restoration of 3D medical images with

total variation scheme on wavelet domains (TVW), in: Proc. SPIE, Vol.

6144, 2006, pp. 61441E–61441E–9.

22

[24] D. C. Adler, T. H. Ko, J. G. Fujimoto, Speckle reduction in optical

coherence tomography images by use of a spatially adaptive wavelet

filter, Opt. Lett. 29 (24) (2004) 2878–2880.

[25] S. Chitchian, M. A. Fiddy, N. M. Fried, Denoising during optical co￾herence tomography of the prostate nerves via wavelet shrinkage us￾ing dual-tree complex wavelet transform, Journal of Biomedical Optics

14 (1) (2009) 014031–014031–6.

[26] Z. Jian, L. Yu, B. Rao, B. J. Tromberg, Z. Chen, Three-dimensional

speckle suppression in optical coherence tomography based on the

curvelet transform, Opt. Express 18 (2) (2010) 1024–1032.

[27] A. Ozcan, A. Bilenca, A. E. Desjardins, B. E. Bouma, G. J. Tearney,

Speckle reduction in optical coherence tomography images using digital

filtering, J. Opt. Soc. Am. A 24 (7) (2007) 1901–1910.

[28] D. L. Marks, T. S. Ralston, S. A. Boppart, Speckle reduction by I￾divergence regularization in optical coherence tomography, J. Opt. Soc.

Am. A 22 (11) (2005) 2366–2371.

[29] A. Wong, A. Mishra, K. Bizheva, D. A. Clausi, General Bayesian es￾timation for speckle noise reduction in optical coherence tomography

retinal imagery, Opt. Express 18 (8) (2010) 8338–8352.

[30] X. Zhang, L. Li, F. Zhu, W. Hou, X. Chen, Spiking cortical model￾based nonlocal means method for speckle reduction in optical coherence

tomography images, Journal of Biomedical Optics 19 (6) (2014) 066005.

23

[31] J. Aum, J. hyun Kim, J. Jeong, Effective speckle noise suppression in

optical coherence tomography images using nonlocal means denoising

filter with double gaussian anisotropic kernels, Appl. Opt. 54 (13) (2015)

D43–D50.

[32] W. Jiang, M. Ding, X. Zhang, Iterative nonlocal means method for

despeckling optical coherence tomography images, Journal of Medical

Imaging and Health Informatics 4 (5) (2014) 819–824.

[33] M. Pircher, E. Gotzinger, R. Leitgeb, A. F. Fercher, C. K. Hitzenberger,

Speckle reduction in optical coherence tomography by frequency com￾pounding, Journal of Biomedical Optics 8 (3) (2003) 565–569.

[34] J. M. Schmitt, S. H. Xiang, K. M. Yung, Speckle in Optical Coherence

Tomography, Journal of Biomedical Optics 4 (1) (1999) 95–105.

[35] A. F. Fercher, W. Drexler, C. K. Hitzenberger, T. Lasser, Optical co￾herence tomography - principles and applications, Reports on Progress

in Physics 66 (2) (2003) 239.

[36] A. Curatolo, B. F. Kennedy, D. D. Sampson, T. Hillman, Speckle in

optical coherence tomography, Taylor & Francis, 2013.

[37] L. Bian, J. Suo, F. Chen, Q. Dai, Multiframe denoising of high-speed op￾tical coherence tomography data using interframe and intraframe priors,

Journal of biomedical optics 20 (3) (2015) 036006–036006.

[38] J. W. Goodman, Some fundamental properties of speckle, J. Opt. Soc.

Am. 66 (11) (1976) 1145–1150.

24

[39] J. Goodman, Statistical Optics, Wiley, New York, 1985.

[40] T. Eltoft, Modeling the amplitude statistics of ultrasonic images, Med￾ical Imaging, IEEE Transactions on 25 (2) (2006) 229–240.

[41] P. M. Shankar, A general statistical model for ultrasonic backscattering

from tissues, Ultrasonics, Ferroelectrics, and Frequency Control, IEEE

Transactions on 47 (3) (2000) 727–736.

[42] M. P. Wachowiak, R. Smol´ıkov´a, J. M. Zurada, A. S. Elmaghraby, Esti￾mation of K distribution parameters using neural networks, IEEE Trans￾actions on Biomedical Engineering 49 (6) (2002) 617–620.

[43] E. E. Kuruo˘glu, J. Zerubia, Modeling SAR images with a generalization

of the Rayleigh distribution, Image Processing, IEEE Transactions on

13 (4) (2004) 527–533.

[44] T. Eltoft, The Rician inverse Gaussian distribution: a new model for

non-Rayleigh signal amplitude statistics, Image Processing, IEEE Trans￾actions on 14 (11) (2005) 1722–1735.

[45] M. P. Shankar, Ultrasonic tissue characterization using a generalized

Nakagami model, Ultrasonics, Ferroelectrics, and Frequency Control,

IEEE Transactions on 48 (6) (2001) 1716–1720.

[46] G. Vegas-S´anchez-Ferrero, D. Mart´ın-Mart´ınez, S. Aja-Fern´andez,

C. Palencia, On the influence of interpolation on probabilistic models

for ultrasonic images, in: Biomedical Imaging: From Nano to Macro,

2010 IEEE International Symposium on, IEEE, 2010, pp. 292–295.

25

[47] C. Tomasi, R. Manduchi, Bilateral filtering for gray and color images,

in: Computer Vision, 1998. Sixth International Conference on, 1998, pp.

839–846.

[48] M. Zhang, B. K. Gunturk, Multiresolution bilateral filtering for image

denoising, Image Processing, IEEE Transactions on 17 (12) (2008) 2324–

2333.

[49] R. Riji, J. Rajan, J. Sijbers, M. S. Nair, Iterative bilateral filter for Ri￾cian noise reduction in MR images, Signal, Image and Video Processing

(2014) 1–6.

[50] I. A. Stegun, M. Abramowitz, Handbook of Mathematical Functions,

with Formulas, Graphs, and Mathematical Tables, Dover Publications,

1964.

[51] A. Clifford Cohen, B. Jones Whitten, Modified moment and maximum

likelihood estimators for parameters of the three-parameter Gamma dis￾tribution, Communications in Statistics-Simulation and Computation

11 (2) (1982) 197–216.

[52] Y. Fiisher, Fractal image compression: Theory and application (1994).

[53] Z. Wang, A. C. Bovik, H. R. Sheikh, E. P. Simoncelli, Image quality as￾sessment: from error visibility to structural similarity, Image Processing,

IEEE Transactions on 13 (4) (2004) 600–612.

[54] A. Bhattachayya, On a measure of divergence between two statistical

population defined by their population distributions, Bulletin Calcutta

Mathematical Society 35 (1943) 99–109.

26

[55] W. K. Pratt (Ed.), Digital Image Processing, Wiley, 1977.

27

Figure 1: Actual distribution of the pixels selected from a homogeneous region of the

image shown in Fig. 5 (a) compared with the Gamma PDF (parameters estimated from

the same homogeneous region)

28

(a)

(b)

Figure 2: Distribution of local estimates of ρ and β from a real swine OCT image (shown

in Fig. 5 (a))(a) Distribution of local estimates of ρ (b) Distribution of local estimates of

β

29

Figure 3: Actual distribution of the pixels selected from a homogeneous region of the

image shown in Fig. 5 (f) compared with the Gaussian PDF (parameters estimated from

the same homogeneous region)

30

(a) (b) (c)

(d) (e) (f) (g)

(h) (i) (j) (k)

(l) (m) (n) (o)

Figure 4: Results obtained with different filters applied to Gamma distributed Phantom image: (a) original Phantom

image; (b)-(c) noisy image (ρ = 4, β = 4) and its residual (d)-(e) Bayesian estimation result and its residual (f)-(g) wavelet

multi-frame filter result and its residual (h)-(i) Wiener filter result and its residual (j)-(k) bilateral filter result and its

residual (l)-(m) proposed AUB filter result and its residual (n)-(o) proposed IAUB filter result and its residual. Residual

images are scaled to 0 − 20.

31

(a)

(b) (c)

(d) (e)

(f) (g)

Figure 5: Results obtained with different filters applied to swine eye OCT image: (a)

original noisy image (b) Bayesian estimation result (c) wavelet multi-frame filter result (d)

Wiener filter result (e) bilateral filter result (f) proposed AUB filter result (g) proposed

IAUB filter result.

32

(a)

(b) (c)

(d) (e)

(f) (g)

Figure 6: Results obtained with different filters applied to Pigs eye OCT image: (a) original

noisy image (b) Bayesian estimation result (c) wavelet multi-frame filter result (d) Wiener

filter result (e) bilateral filter result (f) proposed AUB filter result (g) proposed IAUB

filter result.

33

Table 1: Comparison of results based on PSNR

Filter PSNR(dB)

ρ = 3, β = 2 ρ = 4, β = 3 ρ = 6, β = 6 ρ = 8, β = 7

Bayesian Estimation [29] 21.6847 20.6911 15.7897 12.612

Wavelet Multi-frame [54] 31.3507 25.5902 16.3301 12.6485

Wiener Filtering [16] 31.3993 26.1677 16.9441 13.1653

Bilateral Filtering [47] 32.3345 26.5399 17.0812 13.2393

Proposed AUB Filter 49.5666 45.5476 37.8332 35.5024

Proposed IAUB Filter 53.6759 48.2858 40.6639 38.2314

34

Table 2: Comparison of results based on mean SSIM

Filter mean SSIM

ρ = 3, β = 2 ρ = 4, β = 3 ρ = 6, β = 6 ρ = 8, β = 7

Bayesian Estimation [29] 0.5018 0.4378 0.3545 0.3019

Wavelet Multi-frame [54] 0.5186 0.3908 0.2054 0.1551

Wiener Filtering [16] 0.5581 0.4896 0.4023 0.3479

Bilateral Filtering [47] 0.5721 0.5051 0.4058 0.3261

Proposed AUB Filter 0.9850 0.9692 0.8756 0.8287

Proposed IAUB Filter 0.9920 0.9813 0.9186 0.8830

35

Table 3: Comparison of results based on BC

Filter BC

ρ = 3, β = 2 ρ = 4, β = 3 ρ = 6, β = 6 ρ = 8, β = 7

Bayesian Estimation [29] 0.0347 0.0324 0.0295 0.1155

Wavelet Multi-frame [54] 0.1276 0.0319 0.0766 0.0898

Wiener Filtering [16] 0.0378 0.0326 0.0303 0.1112

Bilateral Filtering [47] 0.0308 0.0145 0.0214 0.1358

Proposed AUB Filter 0.7682 0.7371 0.6291 0.5965

Proposed IAUB Filter 0.8077 0.7877 0.6440 0.5681

36

Table 4: Comparison of results based on Pratt’s FOM

Filter FOM

ρ = 3, β = 2 ρ = 4, β = 3 ρ = 6, β = 6 ρ = 8, β = 7

Bayesian Estimation [29] 0.8404 0.7232 0.4653 0.3590

Wavelet Multi-frame [54] 0.8397 0.7230 0.4651 0.3588

Wiener Filtering [16] 0.8400 0.7240 0.4671 0.3607

Bilateral Filtering [47] 0.8401 0.7243 0.4688 0.3617

Proposed AUB Filter 0.9972 0.9944 0.9859 0.9812

Proposed IAUB Filter 0.9978 0.9954 0.9875 0.9841

37

Table 5: Experts evaluation

Experts Image Type

Grading on filtered result

Bayesian

Estimation

[29]

Wavelet

Multiframe

[54]

Wiener

Filter

[16]

Bilateral

Filter

[47]

Proposed

AUB

Filter

Proposed

IAUB

Filter

Ophthalmologist 1

Swine retinal 3 4 4 4 7 7

Pig eye 3 2 2 3 6 6

Ophthalmologist 2

Swine retinal 6 7 7 5 9 8

Pig eye 6 7 7 6 9 9

38

