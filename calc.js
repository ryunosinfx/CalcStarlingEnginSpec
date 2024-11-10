//
// 初期値設定
//
let lim = 'On';
//
// 物性値の計算
//R:ガス定数(J/kgK)：理想気体の状態方程式 に含まれる気体に特有の定数のことで、単位は [J/(kg・K)] です。
//なお、気体定数に気体の分子量 [kg/mol] をかけたものはすべての理想気体で同一となります。 これを一般気体定数といい、その値は 8.3145 J/(mol・K) となります。
// density[kg/m3]は0℃基準
// 比熱 cp[J/kg ℃]
//
const gasTypes = {
	He: { R: 8.31677 / 0.004003, Myua: 0.0000003709, Myub: 0.6981, density: 0.179, cp: 5192 },
	Air: {
		R: 8.31677 / 0.02895,
		Myua: 0.0000003065, //粘度係数
		Myub: 0.7214,
		density: 1.251,
		cp: 1005,
	},
	N2: {
		R: 8.31677 / 0.028013,
		Myua: 0.0000003326,
		Myub: 0.7005,
		density: 1.211,
		cp: 1043,
	},
	H2: {
		R: 8.31677 / 0.002016,
		Myua: 0.0000001942,
		Myub: 0.6716,
		density: 0.0869,
		cp: 14193,
	},
};
const atom1 = 0.1013;
const M = 1000 * 1000;
const k10 = Math.pow(10, 4);
const ZERO_C = 273;
const Minit = 60;
class Calc {
	static f = null;
	static d = null;
	static init() {
		console.log(document.forms);
		const d = document,
			f = d.forms['in'];
		Calc.d = d;
		Calc.f = f;
		Calc.addInpueEventListener(f.Te);
		Calc.addInpueEventListener(f.Tc);
		Calc.addInpueEventListener(f.Vse);
		Calc.addInpueEventListener(f.Pm);
		for (const gas of f.gas) Calc.addInpueEventListener(gas, 'change');
		Calc.addInpueEventListener(f.cylinderC);
		Calc.addInpueEventListener(f.reV);
		Calc.addInpueEventListener(f.reR);
		Calc.addInpueEventListener(f.reD);
		Calc.addInpueEventListener(f.reCp);
		Calc.addInpueEventListener(f.inT);
		Calc.addInpueEventListener(f.outT);
		Calc.addInpueEventListener(f.leakR);
	}
	static ael = (e, ev, cb, isStopPropagation = false) => (e.addEventListener(ev, cb, isStopPropagation) ? cb : cb);
	static addInpueEventListener(elm, type = 'input') {
		Calc.ael(elm, type, () => Calc.exec());
	}
	/**
	 * 計算値を丸める
	 * @param {*} v 値
	 * @param {*} scale 有効桁数
	 * @param {*} s 開始桁
	 * @param {*} e 終了桁
	 * @returns
	 */
	static rollUp(v, scale = 2, s = -10, e = 0) {
		const sc = Math.pow(10, scale);
		for (let i = s; i <= e; i++) {
			const ii = i + 1,
				k = Math.pow(10, i);
			if (v >= k && v < Math.pow(10, ii)) {
				const v1 = v / k,
					v2 = Math.round(v1 * sc) / sc;
				return v2 + 'E' + i;
			}
		}
	}
	static round = (v, sc = 1) => {
		const p = Math.pow(10, sc);
		return Math.round(v * p) / p;
	};
	static pint = (v) => (!v || isNaN(`${v}`) ? 0 : v * 1);
	static printToId = (id, v) => {
		const elm = Calc.d.getElementById(id);
		elm.textContent = v;
	};

	static exec = (f = Calc.f) => {
		//
		// 計算条件
		//
		const Te = Calc.pint(f.Te.value) + ZERO_C,
			Tc = Calc.pint(f.Tc.value) + ZERO_C,
			V = Calc.pint(f.Vse.value),
			Vse = V / M,
			Pmp = Calc.pint(f.Pm.value),
			cylinderC = Calc.pint(f.cylinderC.value),
			addBC = (Pmp * (1 + Te / Tc)) / 2,
			Pm = Calc.pint(addBC) * M,
			L = Math.pow(Vse, 1 / 3), //容積の一辺
			stroke = Calc.round((L / cylinderC) * 100),
			gas = gasTypes[f.gas.value],
			Myua = gas.Myua,
			Myub = gas.Myub,
			R = gas.R, //ガス定数(J/kgK)
			density = gas.density, // density[kg/m3]は0℃基準
			cp = gas.cp, // 比熱 cp[J/kg ℃]
			vRe = Te / ZERO_C,
			Td = Te - Tc,
			Tv = Te / Tc,
			reV = Calc.pint(f.reV.value), //再生器容積ml
			reR = Calc.pint(f.reR.value) / 100, //再生器：充填率%
			reD = Calc.pint(f.reD.value), //再生器：素材密度g/ml
			reCp = Calc.pint(f.reCp.value), //再生器：素材比熱J/Kg℃
			rff = (reV * reR * reD * reCp) / 1000, //熱容量W※温度差1℃として
			Ete = ((Vse * density) / vRe) * cp * Td * (addBC / atom1),
			Etf = Ete - rff,
			carnotF = 1 - Tc / Te, //カルノー効率
			diameter = Calc.round((L * 100 * 2) / Math.pow(Math.PI, 1 / 2), 4); //直径
		//
		// 無次元量・物理量の計算
		//
		const isLimited = lim === 'On',
			Plim = isLimited ? Pm : Calc.pint(f.Plim.value) * M,
			Tlim = isLimited ? Te : Calc.pint(f.Tlim.value) + ZERO_C,
			vis = ((Myua * Math.pow(Te, Myub)) / Pm) * R * Te, // 動粘性係数
			vislim = ((Myua * Math.pow(Tlim, Myub)) / Plim) * R * Tlim,
			TT = (Te - Tc) / (Te + Tc), //無次元温度T*
			PP = Pm / Plim, //許容圧力比P*
			SS = (Tlim * R * Math.pow(L, 2)) / Math.pow(vislim, 2), //無次元エンジン仕様 S = Tlim*R*Vse^(2/3)/vislim^2
			nn = 0.00023 * Math.pow(SS, 0.56), //無次元回転数nmax=2.3*10^-4*S^0.56
			LLs = 0.15 * Math.pow(nn, 1.06), //最高無次元軸出力Lsmax=0.15*nmax^1.06
			Ns = (nn / Math.pow(L, 2)) * vis, //エンジン回転数
			N = Ns * Minit, //エンジン回転数
			Ls = ((LLs / nn) * Pm * Vse * PP * TT * N) / Minit, //最高軸出力
			Bn = ((Ls / N) * Minit) / Pm / Vse, //ビール数
			Wn = Bn / TT, //ウエスト数
			eFmin = Ls / carnotF,
			outComeF = (eFmin / Calc.pint(f.leakR.value)) * 100,
			inT = Calc.pint(f.inT.value) + ZERO_C,
			outT = Calc.pint(f.outT.value) + ZERO_C,
			Tr = (inT - outT) / outComeF,
			inTr = (inT - Te) / outComeF,
			outTr = (Tc - outT) / eFmin;
		console.log('rff:', rff);
		//
		// 計算結果の表示
		//
		const r = {
			cylinder: Calc.round(L * L * 10000, 4),
			diameter: diameter,
			addVal: Calc.round(V * Tv - V, 4),
			stroke: stroke,
			rodL: Calc.round(stroke * 3, 2),
			minRc: Calc.round(stroke * 3.5, 2),
			cylinderL: Calc.round(stroke + diameter, 2),
			DiffDeg: cylinderC > 1 ? 360 / cylinderC : 90,
			reVpC: Calc.round(reV / cylinderC, 2),
			VsepC: Calc.round(V / cylinderC, 2),
			addBC: Calc.round(addBC, 4),
			carnotF: Calc.round(carnotF * 100, 4),
			bufRq: Calc.round(V * 0.05, 2),
			Tr: Calc.round(Tr, 2),
			inTr: Calc.round(inTr, 2),
			outTr: Calc.round(outTr, 2),
			TeK: Te,
			TcK: Tc,
			visv: Calc.rollUp(vis, 2, -10, 0), // 動粘性係数
			Rv: Calc.round(R), //ガス定数
			SSv: SS < k10 ? Calc.round(SS) : Calc.rollUp(SS, 2, 4, 20), //無次元エンジン仕様
			LLsv: LLs < k10 ? Calc.round(LLs) : Calc.rollUp(LLs, 2, 4, 20), //最高無次元軸出力
			nnv: nn < k10 ? Calc.round(nn) : Calc.rollUp(nn, 2, 4, 20), //無次元回転数
			Lsv: Calc.round(Ls), //最高軸出力
			Nvs: Calc.round(Ns), //エンジン回転数rps
			Nv: Calc.round(N), //エンジン回転数rpm
			eF: Calc.round(Ns * Ete),
			eFF: Calc.round(Ns * Etf),
			eFr: Calc.round((Ls * 100) / (Ns * Etf), 3),
			eFmin: Calc.round(eFmin, 2),
			eFpri: Calc.round(outComeF, 2),
			Bnv: Calc.round(Bn, 3), //ビール数
			Wnv: Calc.round(Wn, 3), //ウエスト数
			rff: Calc.round(rff, 3),
		};
		r.Pc = Calc.round(Ls / ((Ns * L * r.cylinder * M) / 10000), 5);
		r.Pf = Calc.round((r.Pc / (addBC - Pmp)) * 100, 2);
		for (const k in r) Calc.printToId(k, r[k]);
	};
}
Calc.init();
Calc.exec();
