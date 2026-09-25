#include <oxstd.h>
#include <oxprob.h>
#include <packages/ssfpack/ssfpack.h>
//#include <packages/gnudraw/gnudraw.h>
#include <oxdraw.h>
#include "ReportMCMC.ox"
#import <database>

class MFData
{
	MFDiff(const mY, const cN_1);
	MFCenter(const mY, const cN_1);
	MFCenterBreak(const mY, const cN_1, const iBreak);
	MFCalcS2(const mY, const cN_1);
}

MFData::MFDiff(const mY, const cN_1)
{
	decl mdy = diff0(mY[][: cN_1 - 1], 3, .NaN) ~ diff0(mY[][cN_1:], 1, .NaN);
	return mdy;
}

MFData::MFCenter(const mY, const cN_1)
{
	decl vmu = meanc(deleter(mY[][: cN_1 - 1])) ~ meanc(mY[][cN_1:]);
	return mY - vmu;
}

MFData::MFCenterBreak(const mY, const cN_1, const iBreak)
{
	decl my_1 = MFCenter(mY[: iBreak][],    cN_1);
	decl my_2 = MFCenter(mY[iBreak + 1:][], cN_1);
	return my_1 | my_2;
}

MFData::MFCalcS2(const mY, const cN_1)
{
	decl vs2 = (9 / 19) * varc(deleter(mY[][: cN_1 - 1])) ~ varc(mY[][cN_1:]);
	return vs2;
}

class MyMCMC
{
	CalcGamma(const mPhi, const mSigma);
	SSFlarge(const mY, const mPhi, const mSigma, const amA, const amB, const amC, const amD);
	DrawZplus(const mY, const mPhi, const mSigma);
	SSFsmall(const mY, const mPhi, const mSigma, const amLambda, const amA, const amB, const amMu, const amC, const amD);
	CalcInitP(const mY, const mPhi, const mSigma);
	DrawYm(const mY, const mPhi, const mSigma_w);
	DrawSigmaInv(const iN_1, const mS_1Inv);
	DrawPhi(const vM_1, const mSigma, const mD_1Inv);
	CalcLogDens(const vX, const mSigma);
	CalcPmove(const mGamma, const mGamma_star, const vS_0);
	CalcWeights(const mPhi);
}

MyMCMC::CalcGamma(const mPhi, const mSigma)
{
	decl cn     = rows(mPhi);
	decl ip     = columns(mPhi) / cn;
	decl ma     = mPhi | unit((ip - 1) * cn) ~ zeros((ip - 1) * cn, cn);
	decl cm     = rows(ma);
	decl mb     = choleski(mSigma) | zeros(cm - cn, cn);
	decl mgamma = shape(invert(unit(cm ^ 2) - ma ** ma) * vec(mb * mb'), cm, cm);
	return mgamma;
}

MyMCMC::SSFlarge(const mY, const mPhi, const mSigma, const amA, const amB, const amC, const amD)
{
	decl cn   = rows(mY);
	decl cn_1 = rows(selectr(mY));
	decl cn_2 = cn - cn_1;
	decl ip   = columns(mPhi) / cn;

	// Transition equation

	decl ma;
	if (ip <= 4) ma = mPhi ~ zeros(cn, (5 - ip) * cn) | unit(4 * cn) ~ zeros(4 * cn, cn);
	else         ma = mPhi | unit((ip - 1) * cn) ~ zeros((ip - 1) * cn, cn);
	decl cm = rows(ma);
	decl mb = choleski(mSigma) | zeros(cm - cn, cn);

	// Measurement equation

	decl mh_0 = (1 / 3) * unit(cn_1) ~ zeros(cn_1, cn_2) | zeros(cn_2, cn_1) ~ unit(cn_2);
	decl mh_1 = (2 / 3) * unit(cn_1) ~ zeros(cn_1, cn_2) | zeros(cn_2, cn);
	decl mh_2 =           unit(cn_1) ~ zeros(cn_1, cn_2) | zeros(cn_2, cn);
	decl mh_3 = (2 / 3) * unit(cn_1) ~ zeros(cn_1, cn_2) | zeros(cn_2, cn);
	decl mh_4 = (1 / 3) * unit(cn_1) ~ zeros(cn_1, cn_2) | zeros(cn_2, cn);
	decl mc   = mh_0 ~ mh_1 ~ mh_2 ~ mh_3 ~ mh_4;
	if (ip > 5) mc ~= zeros(cn, (ip - 5) * cn);
	decl md = zeros(cn, cn);

	amA[0] = ma;
	amB[0] = mb;
	amC[0] = mc;
	amD[0] = md;
}

MyMCMC::DrawZplus(const mY, const mPhi, const mSigma)
{
	decl ma, mb, mc, md;
	SSFlarge(mY, mPhi, mSigma, &ma, &mb, &mc, &md);
	decl cn     = rows(mY);
	decl ct     = columns(mY);
	decl cm     = rows(ma);
	decl me     = mb | md;
	decl mu     = me * rann(cn, 1 + ct);
	decl mphi   = ma | mc;
	decl momega = me * me';
//	decl mp     = shape(invert(unit(cm ^ 2) - ma ** ma) * vec(mb * mb'), cm, cm);
//	decl mq     = choleski(mp);    // Exact draw
	decl mq     = zeros(cm, cm);   // Approximate draw
	decl mzplus = SsfRecursion(mu, mphi, momega, mq | zeros(1, cm));
	return mzplus;
}

MyMCMC::SSFsmall(const mY, const mPhi, const mSigma, const amLambda, const amA, const amB, const amMu, const amC, const amD)
{
	decl cn      = rows(mY);
	decl cn_1    = rows(selectr(mY));
	decl cn_2    = cn - cn_1;
	decl ip      = columns(mPhi) / cn;
	decl ct      = columns(mY);
	decl mphi_11 = <>;
	decl mphi_12 = <>;
	decl mphi_21 = <>;
	decl mphi_22 = <>;
	decl mx_2    = <>;
	for (decl i = 0; i <= ip - 1; ++i)
	{
		mphi_11 = mphi_11 ~ mPhi[: cn_1 - 1][i * cn : i * cn + cn_1 - 1];
		mphi_12 = mphi_12 ~ mPhi[: cn_1 - 1][i * cn + cn_1 : (i + 1) * cn  - 1];
		mphi_21 = mphi_21 ~ mPhi[cn_1 : cn - 1][i * cn : i * cn + cn_1 - 1];
		mphi_22 = mphi_22 ~ mPhi[cn_1 : cn - 1][i * cn + cn_1 : (i + 1) * cn - 1];
		mx_2    = mx_2 | lag0(mY[cn_1 : cn - 1][]', i + 1)';
	}

	// Transition equation

	decl ma;
	if (ip <= 4) ma = mphi_11 ~ zeros(cn_1, (5 - ip) * cn_1) | unit(4 * cn_1) ~ zeros(4 * cn_1, cn_1);
	else         ma = mphi_11 ~ zeros(cn_1, cn_1) | unit(ip * cn_1) ~ zeros(ip * cn_1, cn_1);
	decl cm      = rows(ma);
	decl mlambda = mphi_12 * mx_2 | zeros(cm - cn_1, ct);
	decl mb      = (unit(cn_1) ~ zeros(cn_1, cn_2)) * choleski(mSigma) | zeros(cm - cn_1, cn);

	// Measurement equation

	decl mmu  = zeros(cn_1, ct) | mphi_22 * mx_2;
	decl mh_0 = (1 / 3) * unit(cn_1);
	decl mh_1 = (2 / 3) * unit(cn_1);
	decl mh_2 = unit(cn_1);
	decl mh_3 = (2 / 3) * unit(cn_1);
	decl mh_4 = (1 / 3) * unit(cn_1);
	decl mh   = mh_0 ~ mh_1 ~ mh_2 ~ mh_3 ~ mh_4;
	decl mc;
	if (ip <= 4) mc = mh | zeros(cn_2, cn_1) ~ mphi_21 ~ zeros(cn_2, (4 - ip) * cn_1);
	else         mc = mh ~ zeros(cn_1, (ip - 4) * cn_1) | zeros(cn_2, cn_1) ~ mphi_21;
	decl md = zeros(cn_1, cn) | (zeros(cn_2, cn_1) ~ unit(cn_2)) * choleski(mSigma);

	// Future form

	decl mmu_star = mmu + mc * mlambda;
	decl mc_star  = mc * ma;
	decl md_star  = mc * mb + md;

	amLambda[0] = mlambda;
	amA[0]      = ma;
	amB[0]      = mb;
	amMu[0]     = mmu_star;
	amC[0]      = mc_star;
	amD[0]      = md_star;
}

MyMCMC::CalcInitP(const mY, const mPhi, const mSigma)
{
	decl cn   = rows(mY);
	decl cn_1 = rows(selectr(mY));
	decl ip   = columns(mPhi) / cn;
	decl ma;
	if (ip <= 4) ma = mPhi ~ zeros(cn, (5 - ip) * cn) | unit(4 * cn) ~ zeros(4 * cn, cn);
	else         ma = mPhi ~ zeros(cn, cn) | unit(ip * cn) ~ zeros(ip * cn, cn);
	decl cm      = rows(ma);
	decl mb      = choleski(mSigma) | zeros(cm - cn, cn);
	decl mp      = shape(invert(unit(cm ^ 2) - ma ** ma) * vec(mb * mb'), cm, cm);
	decl mselect = <>;
	decl mi      = unit(cm);
	for (decl i = 0; i <= cm - 1; ++i)
	{
		if (imod(i, cn) <= cn_1 - 1) mselect |= mi[i][];
	}
	return mselect * mp * mselect';
}

MyMCMC::DrawYm(const mY, const mPhi, const mSigma_w)
{
	decl ma, mb, mc, md;
	SSFlarge(mY, mPhi, mSigma_w, &ma, &mb, &mc, &md);
	decl cn      = rows(mY);
	decl cn_1    = rows(selectr(mY));
	decl ct      = columns(mY);
	decl cm      = rows(ma);
	decl mzplus  = DrawZplus(mY, mPhi, mSigma_w);
	decl myplus  = mzplus[cm : cm + cn - 1][1 :];
	decl mydiff  = mY - myplus;
	decl my1plus = mzplus[: cn_1 - 1][: ct - 1];

	decl mlambda, mmu;
	SSFsmall(mydiff, mPhi, mSigma_w, &mlambda, &ma, &mb, &mmu, &mc, &md);
	cm = rows(ma);
	decl mphi     = ma | mc;
	decl me       = mb | md;
	decl momega   = me * me';
//	decl mp       = CalcInitP(mydiff, mPhi, mSigma_w);              // Unconditional state variance
	decl mp       = 10 ^ 6 * max(1, diagonal(momega)) * unit(cm);   // Diffuse initial state
	decl msigma   = mp | zeros(1, cm);
	decl mdelta   = <>;
	decl mj_phi   = <>;
	decl mj_omega = <>;
	decl mj_delta = range(0, cm + cn - 1)';
	decl mxt      = mlambda | mmu;
	decl mzdiff_s, mudiff_s;
	SsfMomentEst(ST_SMO, &mzdiff_s, mydiff, mphi, momega, msigma, mdelta, mj_phi, mj_omega, mj_delta, mxt);
//	SsfMomentEst(DS_SMO, &mudiff_s, mydiff, mphi, momega, msigma, mdelta, mj_phi, mj_omega, mj_delta, mxt);
	mudiff_s = SsfCondDens(DS_SMO, mydiff, mphi, momega, msigma, mdelta, mj_phi, mj_omega, mj_delta, mxt);
	decl mzdiff_s1 = mxt + mphi * mzdiff_s[: cm  - 1][] + mudiff_s[: cm + cn - 1][];
//	decl mydiff_s  = mxt[cm : cm + cn - 1][] + mzdiff_s[cm : cm + cn - 1][] + mudiff_s[cm : cm + cn - 1][];
//	print(mydiff' ~ mydiff_s');
	decl my1_s     = my1plus + mzdiff_s1[: cn_1 - 1][];
	return my1_s;
}

MyMCMC::DrawSigmaInv(const iN_1, const mS_1Inv)
{
	decl cn        = rows(mS_1Inv);
	decl msigmaInv = choleski(mS_1Inv) * ranwishart(iN_1, cn) * choleski(mS_1Inv)';
	return msigmaInv;
}

MyMCMC::DrawPhi(const vM_1, const mSigma, const mD_1Inv)
{
	decl ck        = rows(vM_1);
	decl cn        = rows(mSigma);
	decl ip        = ck / cn ^ 2;
	decl vmean     = vM_1;
	decl mvariance = mSigma ** mD_1Inv;
	decl vphi, vlambda;
	do
	{
		vphi = vmean + choleski(mvariance) * rann(ck, 1);
		decl mphi = reshape(vphi, cn, ip * cn);
		decl ma   = mphi | unit((ip - 1) * cn) ~ zeros((ip - 1) * cn, cn);
		eigen(ma, &vlambda);
	} while (max(cabs(vlambda)) >= 1);
	return vphi;
}

MyMCMC::CalcLogDens(const vX, const mSigma)
{
	decl dlnf = - (log(determinant(mSigma)) + vX' * invertsym(mSigma) * vX) / 2;
	return dlnf;
}

MyMCMC::CalcPmove(const mGamma, const mGamma_star, const vS_0)
{
	decl dlnf0      = CalcLogDens(vS_0, mGamma);
	decl dlnf0_star = CalcLogDens(vS_0, mGamma_star);
	decl dalpha     = min(exp(dlnf0_star - dlnf0), 1);
	return dalpha;
}

MyMCMC::CalcWeights(const mPhi)
{
	decl cn      = rows(mPhi);
	decl ip      = columns(mPhi) / cn;
	decl ma      = mPhi | unit((ip - 1) * cn) ~ zeros((ip - 1) * cn, cn);
	decl mc_star = unit(cn) ~ zeros(cn, (ip - 1) * cn);
	decl mw      = - mc_star * invert(unit(ip * cn) - ma) * ma;
	return mw;
}

main()
{
	decl time = timer();

	decl dbase = new Database();
	dbase.LoadXls("usq1m3.xls");
	dbase.Info();
	decl mLnX    = dbase.GetAll();
	decl asName  = dbase.GetAllNames();
	decl iYear   = dbase.GetYear1();
	decl iPeriod = dbase.GetPeriod1();
	decl iFreq   = dbase.GetFrequency();
	delete dbase;

	decl cN   = columns(mLnX);
	decl cN_1 = 1;

	decl mfdata = new MFData();
	decl mY = mfdata.MFDiff(mLnX, cN_1)[2:][];
//	mY = mfdata.MFCenter(mY, cN_1);
	decl iBreak = 223;   // Break in 1973M3
	print(mY[iBreak][]);
	mY = mfdata.MFCenterBreak(mY, cN_1, iBreak);   // Comment out if no break
	decl vS2 = mfdata.MFCalcS2(mY, cN_1);
	delete mfdata;

	decl cT = rows(mY);
	decl iP = 12;

	// Prior
	
	decl dLambda = 1;
	decl iN_0    = 6;
	decl mM_0    = zeros(cN, iP * cN);
	decl mD_0    = dLambda ^ -2 * diag(cumulate(ones(iP, 1)) .^ 2) **	diag(vS2);
	decl mD_0Inv = invertsym(mD_0);
	decl mS_0    = (iN_0 - cN - 1) * diag(vS2);

	// Initial value

	decl mphi = mM_0;
	decl msigma = diag(vS2);

	// MCMC

	decl cBurn    = 10000;
	decl cR       = 40000;
	decl mphi_s   = new matrix[cR][iP * cN ^ 2];
	decl msigma_s = new matrix[cR][cN * (cN + 1) / 2];
	decl amy_s    = new array[cN_1];
	decl amgap_s  = new array[cN];
	decl amcorr_s = new array[cN][cN];
	for (decl i = 0; i <= cN - 1; ++i)
	{
		if (i <= cN_1 - 1) amy_s[i] = new matrix[cR][cT];
		amgap_s[i] = new matrix[cR][cT - iP];
		for (decl j = 0; j <= cN - 1; ++j)
		{
			amcorr_s[i][j] = new matrix[cR];
		}
	}
	decl mcmc = new MyMCMC();
//	decl mgamma = mcmc.CalcGamma(mphi, msigma);
	decl mgamma_star;
	for (decl r = -cBurn; r <= cR - 1; ++r)
	{
		// Draw states

		decl my1 = mcmc.DrawYm(mY', mphi, msigma)';
		decl mym = my1 ~ mY[][1 :];

		// OLS

		decl vs_0 = vecr(mym[: iP - 1][]);
		decl my   = mym[iP:][];
		decl ms   = <>;
		decl mx   = <>;
		for (decl i = 1; i <= iP; ++i)
		{
			ms ~= mym[iP + 1 - i: cT - i][];
			mx ~= mym[iP - i: cT - 1 - i][];
		}
		decl mxx      = mx' * mx;
		decl mxxInv   = invertsym(mxx);
		decl mphi_OLS = my' * mx * mxxInv;
		decl me       = my - mx * mphi_OLS';
		decl mee      = me' * me;

		// Draw innovation variance-covariance matrix

		decl in_1           = cT - iP + iN_0;
		decl ms_1           = (mphi_OLS - mM_0) * invertsym(mxxInv + mD_0Inv) * (mphi_OLS - mM_0)' + mee + mS_0;
		decl ms_1Inv        = invertsym(ms_1);
		decl msigmaInv_star = mcmc.DrawSigmaInv(in_1, ms_1Inv);
		decl msigma_star    = invertsym(msigmaInv_star);

		// Draw VAR coefficients

		decl md_1      = mxx + mD_0;
		decl md_1Inv   = invertsym(md_1);
		decl mm_1      = (mphi_OLS * mxx + mM_0 * mD_0) * md_1Inv;
		decl vm_1      = vecr(mm_1);
		decl vphi_star = mcmc.DrawPhi(vm_1, msigma_star, md_1Inv);
		decl mphi_star = reshape(vphi_star, cN, iP * cN);

		// M--H algorithm

//		mgamma_star = mcmc.CalcGamma(mphi_star, msigma_star);
//		decl dalpha = mcmc.CalcPmove(mgamma, mgamma_star, vs_0);
//		if (ranu(1, 1) <= dalpha)
		{
			mphi   = mphi_star;
			msigma = msigma_star;
//			mgamma = mgamma_star;
		}

		// B-N decomposition

		decl mw    = mcmc.CalcWeights(mphi);
		decl mgap  = ms * mw';
		decl mcorr = correlation(mgap);

		// Save results

		if (r >= 0)
		{
			mphi_s[r][]   = vecr(mphi)';
			msigma_s[r][] = vech(msigma)';
			for (decl i = 0; i <= cN - 1; ++i)
			{
				if (i <= cN_1 - 1)  amy_s[i][r][] = mym[][i]';
				amgap_s[i][r][] = mgap[][i]';
				for (decl j = 0; j <= cN - 1; ++j)
				{
					amcorr_s[i][j][r] = mcorr[i][j];
				}
			}
		}
		println(r);
	}

	// Summary

	decl aPhiNames = {};
	for (decl s = 1; s <= iP; ++s)
	{
		for (decl i = 1; i <= cN; ++i)
		{
			for (decl j = 1; j <= cN; ++j)
			{
				aPhiNames ~= sprint("Phi_", s, "(", i, ",", j, ")");
			}
		}
	}
	decl aSigmaNames = {};
	for (decl j = 1; j <= cN; ++j)
	{
		for (decl i = j; i <= cN; ++i)
		{
			aSigmaNames ~= sprint("Sigma(", i, ",", j, ")");
		}
	}

//	decl out = new ReportMCMC(mphi_s);
//	decl out = new ReportMCMC(msigma_s);
	decl out = new ReportMCMC(mphi_s[][: cN ^ 2 - 1]);
	out.SetOutfileName("15f");
	out.SetVarNames(aPhiNames ~ aSigmaNames);
	out.Report();
	delete out;

	// Plots

	decl vinit = mLnX[1][];
	mfdata = new MFData();
	decl my_star = mfdata.MFDiff(mLnX, cN_1)[2:][];
	delete mfdata;
	for (decl i = 0; i <= cN_1 - 1; ++i)
	{
		my_star[][i] = meanc(deleter(my_star[][i])) / 3 + quantilec(amy_s[i])';
		vinit[i]     = mLnX[3][i] - (2 / 3) * my_star[0][i] - (1 / 3) * my_star[1][i];
	}
	delete amy_s;
	decl mar   = vinit + (0 | cumulate(my_star));
	decl amgap = new array[cN];
	decl amnr  = new array[cN];
	decl ampr  = new array[cN];
	for (decl i = 0; i <= cN - 1; ++i)
	{
		amgap[i] = .NaN * ones(3, 1 + iP) ~ quantilec(amgap_s[i], <.5, .05, .95>);
		amnr[i]  = mar[][i]' - amgap[i][1][];
		ampr[i]  = .NaN * ones(1, 1 + iP) ~ meanc(amgap_s[i] .> 0);
	}
	delete amgap_s;

//	SetDraw(SET_COLORMODEL, 1);
	DrawTMatrix(0, mar[][0]' | amnr[0], {"Actual", "Natural"}, iYear, iPeriod, iFreq, 0, 1);
	DrawTMatrix(1, mar[][1]' | amnr[1], {"Actual", "Natural"}, iYear, iPeriod, iFreq, 0, 1);
	DrawTMatrix(2, mar[][2]' | amnr[2], {"Actual", "Natural"}, iYear, iPeriod, iFreq, 0, 1);
	DrawTMatrix(3, mar[][3]' | amnr[3], {"Actual", "Natural"}, iYear, iPeriod, iFreq, 0, 1);
	DrawTitle(0, "Output");
	DrawTitle(1, "Inflation Rate");
	DrawTitle(2, "Nominal Interest Rate");
	DrawTitle(3, "Unemployment Rate");
	DrawAdjust(ADJ_AREAMATRIX, 4, 1);
	SaveDrawWindow("15f_fig4.eps");
//	SaveDrawWindow("15f_sli1.eps");
	ShowDrawWindow();
	DrawTMatrix(0, amgap[0], 0, iYear, iPeriod, iFreq, 0, 1);
	DrawTMatrix(1, amgap[1], 0, iYear, iPeriod, iFreq, 0, 1);
	DrawTMatrix(2, amgap[2], 0, iYear, iPeriod, iFreq, 0, 1);
	DrawTMatrix(3, amgap[3], 0, iYear, iPeriod, iFreq, 0, 1);
	DrawTitle(0, "Output");
	DrawTitle(1, "Inflation Rate");
	DrawTitle(2, "Nominal Interest Rate");
	DrawTitle(3, "Unemployment Rate");
	DrawAdjust(ADJ_AREAMATRIX, 4, 1);
	SaveDrawWindow("15f_fig5.eps");
//	SaveDrawWindow("15f_sli2.eps");
	ShowDrawWindow();
	DrawTMatrix(0, ampr[0], 0, iYear, iPeriod, iFreq, 0, 1);
	DrawTMatrix(1, ampr[1], 0, iYear, iPeriod, iFreq, 0, 1);
	DrawTMatrix(2, ampr[2], 0, iYear, iPeriod, iFreq, 0, 1);
	DrawTMatrix(3, ampr[3], 0, iYear, iPeriod, iFreq, 0, 1);
	DrawTitle(0, "Output");
	DrawTitle(1, "Inflation Rate");
	DrawTitle(2, "Nominal Interest Rate");
	DrawTitle(3, "Unemployment Rate");
	DrawAdjust(ADJ_AREAMATRIX, 4, 1);
	SaveDrawWindow("15f_fig6.eps");
//	SaveDrawWindow("15f_sli3.eps");
	ShowDrawWindow();
	DrawXMatrix(0, amgap[1][0][], "", amgap[0][0][], "", 1, 1);
	DrawXMatrix(2, amgap[1][0][], "", amgap[3][0][], "", 1, 1);
	DrawXMatrix(4, amgap[0][0][], "", amgap[3][0][], "", 1, 1);
	DrawText(0, "Output Gap",            0, 0, -1, -1, TEXT_XLABEL);
	DrawText(0, "Inflation Rate Gap",    0, 0, -1, -1, TEXT_YLABEL);
	DrawText(2, "Unemployment Rate Gap", 0, 0, -1, -1, TEXT_XLABEL);
	DrawText(2, "Inflation Rate Gap",    0, 0, -1, -1, TEXT_YLABEL);
	DrawText(4, "Unemployment Rate Gap", 0, 0, -1, -1, TEXT_XLABEL);
	DrawText(4, "Output Gap",            0, 0, -1, -1, TEXT_YLABEL);
	DrawDensity(1, amcorr_s[0][1][]', "", 1, 0, 0, 0, 0, 0, 1);
	DrawDensity(3, amcorr_s[3][1][]', "", 1, 0, 0, 0, 0, 0, 1);
	DrawDensity(5, amcorr_s[0][3][]', "", 1, 0, 0, 0, 0, 0, 1);
	DrawAdjust(ADJ_AREA_X, 1, -1, 1);
	DrawAdjust(ADJ_AREA_X, 3, -1, 1);
	DrawAdjust(ADJ_AREA_X, 5, -1, 1);
	DrawTitle(1, "Correlation Coefficient");
	DrawTitle(3, "Correlation Coefficient");
	DrawTitle(5, "Correlation Coefficient");
	DrawAdjust(ADJ_AREAMATRIX, 3, 2);
	SaveDrawWindow("15f_fig7.eps");
	ShowDrawWindow();
	DrawDensity(0, amcorr_s[0][1][]', "", 1, 0, 0, 0, 0, 0, 1);
	DrawAdjust(ADJ_AREA_X, 0, -1, 1);
	DrawTitle(0, "Correlation Coefficient");
//	SaveDrawWindow("15f_sli4.eps");
	ShowDrawWindow();
	DrawDensity(0, amcorr_s[0][3][]', "", 1, 0, 0, 0, 0, 0, 1);
	DrawAdjust(ADJ_AREA_X, 0, -1, 1);
	DrawTitle(0, "Correlation Coefficient");
//	SaveDrawWindow("15f_sli5.eps");
	ShowDrawWindow();

	decl vmrgdp   = mar[][0];
	decl vmgdpgap = amgap[0][0][]';
	print(vmrgdp ~ vmgdpgap);
	savemat("temp.xls", vmrgdp ~ vmgdpgap);

	println(timespan(time, timer()));
}
