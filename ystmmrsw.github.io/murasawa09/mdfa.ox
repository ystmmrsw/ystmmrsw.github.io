// This program outputs the followings:
// (i) EMD or OMD estimates of the parameters with asymptotic standard errors,
// (ii) distance test statistic (valid with OMD estimates) and the p-value,
// (iii) weights for the factor scores with asymptotic standard errors.
//
// Caution: This program does not impose implicit parameter restrictions
// (positive definiteness of variance-covariance matrices, covariance stationarity, etc.).
// Make sure that the estimates satisfy such restrictions.

#include <oxstd.h>
#import <maximize>

class Sample
{
	autocov(const mX, const iS, const iLag);
	gammahat(const mX, const iS);
}

Sample::autocov(const mX, const iS, const iLag)
{
	decl ct = rows(mX);
	decl mx = mX[iS: ct - 1][];
	decl mxlag = mX[iS - iLag: ct - 1 - iLag][];
	return (mx - meanc(mx))' * (mxlag - meanc(mx)) / (ct - iS);
}

Sample::gammahat(const mX, const iS)
{
	decl mgamma = autocov(mX, iS, 0);
	decl vgamma = vech(mgamma);
	if (iS > 0)
	{
		decl s;
		for (s = 1; s <= iS; ++s)
		{
			mgamma = autocov(mX, iS, s);
			vgamma = vgamma | vec(mgamma);
		}
	}
	return vgamma;
}

class HAC
{
	x2z(const mX, const iS);
	prewhiten(const mX, const amPhi, const amE);
	autocov(const mX, const iLag);
	NW(const mX);
	recolor(const mSigma, const mPhi);
}

HAC::x2z(const mX, const iS)
{
	decl ct = rows(mX);
	decl mx = mX - meanc(mX[iS: ct - 1][]);
	mx = mx';
	decl mz = vech(mx[][iS] * mx[][iS]');
	if (iS > 0)
	{
		decl s;
		for (s = 1; s <= iS; ++s)
		{
			mz = mz | vec(mx[][iS] * mx[][iS - s]');
		}
	}
	decl vz;
	decl t;
	for (t = iS + 1; t < ct; ++t)
	{
		vz = vech(mx[][t] * mx[][t]');
		if (iS > 0)
		{
			decl s;
			for (s = 1; s <= iS; ++s)
			{
				vz = vz | vec(mx[][t] * mx[][t - s]');
			}
		}
	mz = mz ~ vz;
	}
	return mz';
}

HAC::prewhiten(const mX, const amPhi, const amE)
{
	decl ct = rows(mX);
	decl mx = mX - meanc(mX);
	decl mxlag = lag0(mx, 1);
	decl mphi = mx' * mxlag * invertsym(mxlag' * mxlag);
	decl me = mx - mxlag * mphi';
	amPhi[0] = mphi;
	amE[0] = me[1: ct - 1][];
}

HAC::autocov(const mX, const iLag)
{
	decl ct = rows(mX);
	decl mx = mX - meanc(mX);
	decl mxlag = lag0(mx, iLag);
	return mx' * mxlag / ct;
}

HAC::NW(const mX)
{
	decl ct = rows(mX);
	decl cn = columns(mX);
	decl vw = ones(cn, 1);
	decl vx = mX * vw;
	decl dn = 4 * (ct / 100) ^ (2 / 9);
	decl ds0 = varc(vx);
	decl ds1 = 0;
	decl s;
	for (s = 1; s < dn; ++s)
	{
		ds0 = ds0 + 2 * autocov(vx, s);
		ds1 = ds1 + s * autocov(vx, s);
	}
	decl dgamma = 1.1447 * ((ds1 / ds0) ^ 2) ^ (1 / 3);
	decl dm = dgamma * ct ^ (1 / 3);
	decl msigma = variance(mX);
	for (s = 1; s < dm; ++s)
	{
		msigma = msigma + (1 - s / ceil(dm)) * (autocov(mX, s) + autocov(mX, s)');
	}
	return msigma;
}

HAC::recolor(const mSigma, const mPhi)
{
	decl cn = rows(mSigma);
	return invert(unit(cn) - mPhi) * mSigma * invert(unit(cn) - mPhi)';
}

class MDFA
{
	decl m_cT, m_cN, m_cK, m_iS, m_vGammahat, m_mW;
	MDFA(const cT, const cN, const cK, const iS, const vGammahat, const mW);
	guess(const cN, const cK, const iS);
	theta2gamma(const avGamma, const vTheta);
	distance(const vTheta, const adFunc, const avScore, const amHessian);
	estimate(const cN, const cK, const iS, const avTheta, const adTest);
	asycov(const cT, const iS, const vTheta, const mW, const mSigma);
}

MDFA::MDFA(const cT, const cN, const cK, const iS, const vGammahat, const mW)
{
	m_cT = cT;
	m_cN = cN;
	m_cK = cK;
	m_iS = iS;
	m_vGammahat = vGammahat;
	m_mW = mW;
}

MDFA::guess(const cN, const cK, const iS)
{
	decl mb_2 = ones(cN - cK, cK);
	decl mgamma_ff = unit(cK);
	decl vgamma_uu = ones(cN, 1);
	decl vtheta = vec(mb_2) | vech(mgamma_ff) | vgamma_uu;
	if (iS > 0)
	{
		decl s;
		for (s = 1; s <= iS; ++s)
		{
			mgamma_ff = zeros(cK, cK);
			vgamma_uu = zeros(cN, 1);
			vtheta = vtheta | vec(mgamma_ff) | vgamma_uu;
		}
	}
	return vtheta;
}

MDFA::theta2gamma(const avGamma, const vTheta)
{
	decl cn = m_cN;
	decl ck = m_cK;
	decl is = m_iS;
	decl cb_2 = (cn - ck) * ck;
	decl cgamma_ff = ck * (ck + 1) / 2;
	decl cgamma_uu = cn;
	decl cpara = cb_2 + cgamma_ff + cgamma_uu;
	decl mb_2 = shape(vTheta[0: cb_2 - 1], cn - ck, ck);
	decl mb = unit(ck) | mb_2;
	decl mgamma_ff = unvech(vTheta[cb_2: cb_2 + cgamma_ff - 1]);
	decl mgamma_uu = diag(vTheta[cb_2 + cgamma_ff: cpara - 1]);
	decl mgamma_xx = mb * mgamma_ff * mb' + mgamma_uu;
	decl vgamma_xx = vech(mgamma_xx);
	if (is > 0)
	{
		decl s;
		for (s = 1; s <= is; ++s)
		{
			mgamma_ff = shape(vTheta[cpara: cpara + ck ^ 2 - 1], ck, ck);
			mgamma_uu = diag(vTheta[cpara + ck ^ 2: cpara + ck ^ 2 + cn - 1]);
			mgamma_xx = mb * mgamma_ff * mb' + mgamma_uu;
			vgamma_xx = vgamma_xx | vec(mgamma_xx);
			cpara = cpara + ck ^ 2 + cn;
		}
	}
	avGamma[0] = vgamma_xx;
	return 1;
}

MDFA::distance(const vTheta, const adFunc, const avScore, const amHessian)
{
	decl ct = m_cT;
	decl is = m_iS;
	decl vgammahat = m_vGammahat;
	decl mw = m_mW;
	decl vgamma;
	theta2gamma(&vgamma, vTheta);
	adFunc[0] = - (ct - is) * (vgammahat - vgamma)' * mw * (vgammahat - vgamma);
	return 1;
}

MDFA::estimate(const cN, const cK, const iS, const avTheta, const adTest)
{
	decl vtheta = guess(cN, cK, iS);
	decl dfunc;
	MaxControl(-1, -1);
	MaxBFGS(distance, &vtheta, &dfunc, 0, 1);
	avTheta[0] = vtheta;
	adTest[0] = - dfunc;
}

MDFA::asycov(const cT, const iS, const vTheta, const mW, const mSigma)
{
	decl mj;
	NumJacobian(theta2gamma, vTheta, &mj);
	return invertsym(mj' * mW * mj) * mj' * mW * mSigma * mW * mj * invertsym(mj' * mW * mj) / (cT - iS);
}

class Weights
{
	decl m_cN, m_cK;
	Weights(const cN, const cK);
	theta2weights(const amWeights, const vTheta);
	asycov(const vTheta, const mAsycov);
}

Weights::Weights(const cN, const cK)
{
	m_cN = cN;
	m_cK = cK;
}

Weights::theta2weights(const avWeights, const vTheta)
{
	decl cn = m_cN;
	decl ck = m_cK;
	decl cb_2 = (cn - ck) * ck;
	decl cgamma_ff = ck * (ck + 1) / 2;
	decl cgamma_uu = cn;
	decl cpara = cb_2 + cgamma_ff + cgamma_uu;
	decl mb_2 = shape(vTheta[0: cb_2 - 1], cn - ck, ck);
	decl mb = unit(ck) | mb_2;
	decl mgamma_uu = diag(vTheta[cb_2 + cgamma_ff: cpara - 1]);
	decl mweights = invertsym(mb' * invertsym(mgamma_uu) * mb) * mb' * invertsym(mgamma_uu);
	mweights = mweights ./ sumr(mweights);   // Make weights add up to 1 if necessary.
	avWeights[0] = vec(mweights');
	return 1;
}

Weights::asycov(const vTheta, const mAsycov)
{
	decl mj;
	NumJacobian(theta2weights, vTheta, &mj);
	return mj * mAsycov * mj';
}

main()
{
	decl mX = loadmat("bci.xls");
	mX = standardize(mX);   // Standardize the series if necessary.
	decl ct = rows(mX);
	decl cn = columns(mX);
	decl cK = 1;
	decl iS = 0;
	decl cmome = cn * (cn + 1) / 2 + iS * cn ^ 2;
	decl cpara = (cn - cK) * cK + cK * (cK + 1) / 2 + cn + iS * (cK ^ 2 + cn);

	decl sampleobj = new Sample();
	decl vgammahat = sampleobj.gammahat(mX, iS);
	delete sampleobj;
	
	decl hacobj = new HAC();
	decl mz = hacobj.x2z(mX, iS);	
	decl mphi, me;
	hacobj.prewhiten(mz, &mphi, &me);
	decl msigma = hacobj.NW(me);
	msigma = hacobj.recolor(msigma, mphi);
	delete hacobj;

	decl mw = unit(rows(vgammahat));
	mw = invertsym(msigma);   // Only for the OMD estimator.
	
	decl mdfaobj = new MDFA(ct, cn, cK, iS, vgammahat, mw);
	decl vtheta, dtest;
	mdfaobj.estimate(cn, cK, iS, &vtheta, &dtest);
	decl masycov = mdfaobj.asycov(ct, iS, vtheta, mw, msigma);
	delete mdfaobj;

	decl weightsobj = new Weights(cn, cK);
	decl mweights;
	weightsobj.theta2weights(&mweights, vtheta);
	decl masycovw =weightsobj.asycov(vtheta, masycov);
	delete weightsobj;
	
	print(vgammahat, msigma);
	print(vtheta, sqrt(diagonal(masycov))');
	print(dtest, cmome - cpara, 1 - probchi(dtest, cmome - cpara));
	print(mweights, sqrt(diagonal(masycovw))');
}
