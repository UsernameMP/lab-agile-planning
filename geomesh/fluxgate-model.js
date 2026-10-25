/* GeoMesh Fluxgate v2: single-axis Vacquier R&D model.
   Plain JS for GitHub Pages and Node tests. Nothing here is a bench-validated
   materials/noise/drift model. Benchmarked results are NOT fitted parameters.
*/
(function(root,factory){
  const api=factory();
  if(typeof module==='object'&&module.exports)module.exports=api;
  else root.GeoMeshFluxgate=api;
})(typeof globalThis!=='undefined'?globalThis:this,function(){
  'use strict';
  const MU0=4*Math.PI*1e-7, K_B=1.380649e-23, RHO_CU=1.724e-8, T0=293.15;
  const B_TEST=50e-6, B_FULL=100e-6, TEMPERATURE_RESIDUAL_K=2, GAIN_DRIFT_PPM_PER_K=2;
  // All property values below are design hypotheses / material-class priors.
  // In particular, nRef, age24, temp are NOT verified catalog data.
  const MATERIALS={
    ni:{label:'Ni80Fe20 · 50 мкм, 2 слоя',kind:'foil',tMm:.050,layers:2,mu:38000,hk:160,
      nRef:36,age24:.70,temp:.26,referenceAreaMm2:.20,referenceLengthMm:45,referenceFreqKhz:12,
      width:{min:1,max:3.5,step:.1,def:2.0},status:'гипотеза, параметры партии не измерены'},
    mo:{label:'NiFeMo / 1J85 · 25 мкм, 4 слоя',kind:'foil',tMm:.025,layers:4,mu:68000,hk:120,
      nRef:29,age24:.55,temp:.20,referenceAreaMm2:.20,referenceLengthMm:45,referenceFreqKhz:12,
      width:{min:1,max:3.5,step:.1,def:2.0},status:'гипотеза, отжиг и состав не подтверждены'},
    co:{label:'CoFeSiB · аморфная проволока',kind:'wire',tMm:null,layers:1,mu:27000,hk:120,
      nRef:26,age24:1.15,temp:.50,referenceAreaMm2:Math.PI*.16**2/4,referenceLengthMm:45,referenceFreqKhz:12,
      width:{min:.10,max:.25,step:.01,def:.16},status:'параллельная проволочная схема; иные статьи не переносить напрямую'}
  };
  const FIXED={
    turnExcDensityPerMm:250/45,excWireMm:.22,
    compensatorMaxA:.020,compensatorWireMm:.15,
    coilFillPitch:1.12,maxPickupLayers:8,
    analogueMw:8,adcMw:5,clockMw:2,driverEfficiency:.65,
    inputVoltageNoiseNv:3.0,             // white input voltage noise: nV/sqrt(Hz)
    tempResidualK:TEMPERATURE_RESIDUAL_K,
    gainDriftPpmPerK:GAIN_DRIFT_PPM_PER_K,
    compensationTestT:B_TEST,fullscaleT:B_FULL,
    // All environmental errors are globally fixed; tilt & XYZ do not enter head optimization.
    maxLengthMm:70,maxOuterDiameterMm:35,maxDriveVoltageV:3.3,maxPowerMw:160
  };
  const DEFAULT={material:'ni',lengthMm:45,widthMm:2,spacingMm:6,pickupTurns:1600,freqKhz:12,driveRatio:2.7};
  const LIMITS={lengthMm:[30,70],spacingMm:[4,12],pickupTurns:[600,2400],freqKhz:[6,28],driveRatio:[1.5,4.5]};
  const clamp=(v,min,max)=>Math.max(min,Math.min(max,v));
  const stable=(v)=>Number.isFinite(v)&&v>=0;
  function canonical(src){
    const x=Object.assign({},DEFAULT,src||{});
    x.material=MATERIALS[x.material]?x.material:DEFAULT.material;
    x.widthMm=clamp(+x.widthMm||MATERIALS[x.material].width.def,MATERIALS[x.material].width.min,MATERIALS[x.material].width.max);
    for(const [key,[lo,hi]] of Object.entries(LIMITS))x[key]=clamp(+x[key]||DEFAULT[key],lo,hi);
    x.pickupTurns=50*Math.round(x.pickupTurns/50);
    return x;
  }
  function calc(params){
    const p=canonical(params),m=MATERIALS[p.material],rad=Math.PI/180;
    const L=p.lengthMm*1e-3, A=m.kind==='foil'?p.widthMm*m.tMm*m.layers*1e-6:Math.PI*Math.pow(p.widthMm*1e-3/2,2);
    const areaMm2=A*1e6, V=A*L, referenceV=m.referenceAreaMm2*1e-6*m.referenceLengthMm*1e-3;
    // Proxy for an open-ended thin rod: equivalent-radius prolate shape; not a 3D-FEM result.
    const re=Math.sqrt(A/Math.PI);
    const ratio=Math.max(2.1,2*L/re);
    const D=clamp(Math.pow(re/L,2)*(Math.log(ratio)-1),1e-6,.3);
    const muApp=(H)=>{
      const mr=1+(m.mu-1)/(1+Math.pow(Math.abs(H/m.hk),4));
      return mr/(1+D*(mr-1));
    };
    const excitationTurns=clamp(Math.round(p.lengthMm*FIXED.turnExcDensityPerMm/10)*10,160,400);
    const Hpk=p.driveRatio*m.hk, Ipeak=Hpk*L/excitationTurns;
    const N=128; let mu2=0;
    const wave=[];
    for(let i=0;i<N;i++){
      const phi=2*Math.PI*i/N,v=muApp(Hpk*Math.sin(phi));
      mu2+=v*Math.cos(2*phi);wave.push(v);
    }
    mu2=2*mu2/N;
    const f=p.freqKhz*1000;
    // External field in Tesla -> second-harmonic induced voltage in V.
    // Two opposite-driven rods add external-field modulation in a common pickup.
    const sensitivityVperT=4*Math.PI*f*(2*p.pickupTurns*A)*Math.abs(mu2);
    const copperOhmPerMeter=(diamMm)=>RHO_CU/(Math.PI*Math.pow(diamMm*1e-3/2,2));
    const pickupWireStock=[.20,.18,.16,.15,.14,.12,.10,.08];
    let pickupWireMm=.08,pickupLayers=0,axialTurns=0;
    for(const dia of pickupWireStock){
      const perLayer=Math.max(1,Math.floor(.85*p.lengthMm/(FIXED.coilFillPitch*dia)));
      const layers=Math.ceil(p.pickupTurns/perLayer);
      if(layers<=FIXED.maxPickupLayers){pickupWireMm=dia;pickupLayers=layers;axialTurns=perLayer;break;}
    }
    if(!pickupLayers){axialTurns=Math.floor(.85*p.lengthMm/(FIXED.coilFillPitch*pickupWireMm));pickupLayers=Math.ceil(p.pickupTurns/axialTurns);}
    const halfWidth=m.kind==='wire'?p.widthMm/2:p.widthMm/2;
    const coreRadiusMm=halfWidth+.48;
    const pickupRadiusMm=(p.spacingMm/2)+coreRadiusMm+pickupLayers*pickupWireMm*.56+.25;
    const pickupR=copperOhmPerMeter(pickupWireMm)*(2*Math.PI*pickupRadiusMm*1e-3)*p.pickupTurns*1.12;
    const excitationRadiusMm=Math.sqrt(areaMm2/Math.PI)+.65+FIXED.excWireMm*.5;
    const excitationR=copperOhmPerMeter(FIXED.excWireMm)*2*Math.PI*excitationRadiusMm*1e-3*excitationTurns*1.15;
    const excitationCopperMw=1e3*2*(Ipeak*Ipeak/2)*excitationR;
    // Closed-loop feedback: separate long winding around the pair.
    const feedbackLengthM=L*1.24;
    const targetTurns=B_FULL*feedbackLengthM/(MU0*FIXED.compensatorMaxA);
    const feedbackTurns=clamp(25*Math.ceil(targetTurns/25),200,1000);
    const feedbackIatBtest=B_TEST*feedbackLengthM/(MU0*feedbackTurns);
    const feedbackR=copperOhmPerMeter(FIXED.compensatorWireMm)*2*Math.PI*(pickupRadiusMm+1.0)*1e-3*feedbackTurns;
    const feedbackCopperMw=1e3*feedbackIatBtest**2*feedbackR;
    // Prior core loss is only a parametric engineering allowance: no B-H loop supplied.
    const coreLossMw=2*1.6*Math.pow(Math.max(.1,V/referenceV),.65)*(p.freqKhz/12)*Math.pow(p.driveRatio/2.7,1.6);
    const powerMw=FIXED.analogueMw+FIXED.adcMw+FIXED.clockMw+(excitationCopperMw+coreLossMw)/FIXED.driverEfficiency+feedbackCopperMw;
    const Uexc=Ipeak*Math.hypot(excitationR,2*Math.PI*f*(MU0*60*excitationTurns**2*A/L));
    const outerDiameterMm=2*(pickupRadiusMm+1.6);
    const // Phenomenological fluxgate core-noise prior; NOT derivable from geometry without measurements.
      volumeFactor=Math.pow(referenceV/V,.28),
      freqFactor=Math.pow(12/p.freqKhz,.22),
      drivePenalty=1+.35*Math.pow(p.driveRatio/2.7-1,2),
      spacingPenalty=1+.07*Math.pow(Math.max(0,6-p.spacingMm)/2,2),
      highFreqPenalty=1+.15*Math.pow(Math.max(0,p.freqKhz-18)/10,2);
    const magneticAsd1=m.nRef*volumeFactor*freqFactor*drivePenalty*spacingPenalty*highFreqPenalty;
    const inputEn=Math.hypot(Math.sqrt(4*K_B*T0*pickupR),FIXED.inputVoltageNoiseNv*1e-9);
    const electronicAsd1=inputEn/Math.max(sensitivityVperT,1e-12)*1e12; // pT/sqrtHz
    const asd1=Math.hypot(magneticAsd1,electronicAsd1);
    // 24-hour drift proxy: material-class 1-sigma aging, residual 2K thermal variation,
    // feedback gain drift on 50 uT. No assertion of measured Allan deviation.
    const aging24=m.age24*Math.pow(referenceV/V,.18)*(1+.12*(p.lengthMm/45-1)**2);
    const thermal24=m.temp*FIXED.tempResidualK*
      (1+.09*Math.abs(p.driveRatio/2.7-1)+.06*Math.max(0,(6-p.spacingMm)/2));
    const gain24=(FIXED.compensationTestT*1e9)*FIXED.gainDriftPpmPerK*1e-6*FIXED.tempResidualK;
    const drift24=Math.hypot(aging24,thermal24,gain24);
    const pass=(
      p.driveRatio>=1.6 && p.driveRatio<=4.4 && Ipeak<=.17 && excitationTurns<=400 &&
      pickupLayers<=FIXED.maxPickupLayers && outerDiameterMm<=FIXED.maxOuterDiameterMm &&
      p.lengthMm<=FIXED.maxLengthMm && Uexc<=FIXED.maxDriveVoltageV &&
      feedbackIatBtest*2<=FIXED.compensatorMaxA && powerMw<=FIXED.maxPowerMw &&
      Number.isFinite(asd1)&&Number.isFinite(drift24)
    );
    const flags=[
      ['excitation',Ipeak<=.17&&Uexc<=FIXED.maxDriveVoltageV],
      ['winding',pickupLayers<=FIXED.maxPickupLayers],
      ['envelope',outerDiameterMm<=FIXED.maxOuterDiameterMm],
      ['compensation',feedbackIatBtest*2<=FIXED.compensatorMaxA],
      ['saturation',p.driveRatio>=1.6&&p.driveRatio<=4.4]
    ].filter(x=>!x[1]).map(x=>x[0]);
    return {input:p,material:m,volumeMm3:V*1e9,areaMm2,demagD:D,mu2:Math.abs(mu2),muAppZero:muApp(0),
      Hpk,Ipeak,excitationTurns,excitationR,excitationCopperMw,Uexc,
      pickupWireMm,pickupR,pickupLayers,pickupAxialTurns:axialTurns,pickupRadiusMm,
      feedbackTurns,feedbackR,feedbackIatBtest,feedbackCopperMw,
      sensitivityVperT,magneticAsd1,electronicAsd1,asd1,
      aging24,thermal24,gain24,drift24,coreLossMw,powerMw,outerDiameterMm,pass,flags,wave};
  }
  function dominates(a,b){
    // weakly no worse in each criterion, strictly better in at least one.
    return a.asd1<=b.asd1&&a.drift24<=b.drift24&&a.powerMw<=b.powerMw&&
      (a.asd1<b.asd1||a.drift24<b.drift24||a.powerMw<b.powerMw);
  }
  function pareto(list){
    const order=list.slice().sort((a,b)=>a.asd1-b.asd1),front=[];
    for(const x of order){
      if(front.some(y=>dominates(y,x)))continue;
      for(let i=front.length-1;i>=0;i--)if(dominates(x,front[i]))front.splice(i,1);
      front.push(x);
    }
    return front.sort((a,b)=>a.powerMw-b.powerMw);
  }
  function explore(material,count=5000,seed=9172026){
    const m=MATERIALS[material]||MATERIALS.ni;
    let state=seed>>>0;
    function rand(){state=(1664525*state+1013904223)>>>0;return state/4294967296;}
    const lerp=(a,b)=>a+(b-a)*rand();
    const results=[],rejected={};
    for(let i=0;i<count;i++){
      const p={material:material in MATERIALS?material:'ni',
        lengthMm:Math.round(lerp(30,70)),
        widthMm:Math.round(lerp(m.width.min,m.width.max)/m.width.step)*m.width.step,
        spacingMm:Math.round(lerp(4,12)*2)/2,
        pickupTurns:Math.round(lerp(600,2400)/50)*50,
        freqKhz:Math.round(lerp(6,28)),
        driveRatio:Math.round(lerp(1.6,4.4)*10)/10};
      const r=calc(p);
      if(r.pass)results.push(r);
      else for(const f of r.flags)rejected[f]=(rejected[f]||0)+1;
    }
    const front=pareto(results);
    return {results,front,attempted:count,valid:results.length,rejected,seed,material};
  }
  return {MATERIALS,FIXED,DEFAULT,LIMITS,canonical,calc,dominates,pareto,explore};
});
