'use strict';
const test=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const M=require('../fluxgate-model.js');
const near=(a,b,eps=1e-8)=>Math.abs(a-b)<=eps*Math.max(1,Math.abs(a),Math.abs(b));

test('nominal one-axis dual-core design passes physical constraints',()=>{
 const r=M.calc(M.DEFAULT);
 assert.equal(r.pass,true,JSON.stringify(r.flags));
 assert.ok(r.Ipeak>0&&r.Ipeak<.17);
 assert.ok(r.asd1>0&&r.powerMw>0&&r.drift24>0);
 assert.ok(r.sensitivityVperT>0);
 assert.ok(near(r.Ipeak*r.excitationTurns/(r.input.lengthMm/1000),r.input.driveRatio*r.material.hk));
 assert.ok(near(r.areaMm2,.2));
 assert.ok(r.feedbackIatBtest*2<=M.FIXED.compensatorMaxA);
 assert.ok(r.pickupLayers<=M.FIXED.maxPickupLayers);
});
test('parameter categories are distinct: independent design vs derived physics',()=>{
 assert.equal(Object.keys(M.DEFAULT).length,7); // 1 material, 6 design variables
 const r=M.calc(M.DEFAULT);
 assert.ok(r.pickupWireMm>0 && r.feedbackTurns>=200);
 assert.ok(r.excitationTurns>=160);
 assert.ok(!Object.hasOwn(M.DEFAULT,'gain24'));
 assert.ok(!Object.hasOwn(M.DEFAULT,'vertical'));
 assert.ok(!Object.hasOwn(M.DEFAULT,'tilt'));
});
test('noise is a hypothesis, not borrowed from 5pT benchmark',()=>{
 const r=M.calc(M.DEFAULT);
 assert.ok(near(r.material.nRef,36));
 assert.ok(r.asd1>=r.magneticAsd1);
 assert.ok(r.electronicAsd1>=0);
 assert.ok(r.material.status.includes('гипотеза'));
});
test('reference material alternatives produce distinct physical areas',()=>{
 const ni=M.calc(M.DEFAULT);
 const co=M.calc({...M.DEFAULT,material:'co',widthMm:.16});
 assert.equal(co.material.kind,'wire');
 assert.ok(co.areaMm2<ni.areaMm2);
 assert.ok(ni.input.lengthMm===co.input.lengthMm);
 assert.ok(co.asd1>0&&co.drift24>0);
});
test('current follows saturation reserve and power accounts for copper',()=>{
 const a=M.calc({...M.DEFAULT,driveRatio:2});
 const b=M.calc({...M.DEFAULT,driveRatio:3});
 assert.ok(b.Ipeak>a.Ipeak);
 assert.ok(b.excitationCopperMw>a.excitationCopperMw);
 assert.ok(b.powerMw>a.powerMw);
});
test('pickup turns change electrical slope (not independent measured noise)',()=>{
 const a=M.calc({...M.DEFAULT,pickupTurns:1000});
 const b=M.calc({...M.DEFAULT,pickupTurns:1900});
 assert.ok(b.sensitivityVperT>a.sensitivityVperT);
 assert.ok(b.pickupR>a.pickupR);
});
test('24-hour zero offset drift excludes feedback gain drift',()=>{
 const r=M.calc(M.DEFAULT);
 assert.ok(near(r.drift24,Math.hypot(r.aging24,r.thermal24)));
 assert.ok(near(r.readingDrift24,Math.hypot(r.drift24,r.gain24)));
 assert.ok(r.readingDrift24>r.drift24);
});
test('dominance compares all three objectives',()=>{
 const a={asd1:10,drift24:5,powerMw:30};
 const b={asd1:11,drift24:6,powerMw:40};
 const c={asd1:9,drift24:9,powerMw:20};
 assert.equal(M.dominates(a,b),true);
 assert.equal(M.dominates(a,c),false);
 assert.equal(M.dominates(c,a),false);
 assert.equal(M.pareto([a,b,c]).length,2);
});
test('5000 candidates are generated deterministically and Pareto is truly nondominated',()=>{
 const run=M.explore('ni',5000);
 assert.equal(run.attempted,5000);
 assert.ok(run.valid>100);
 assert.ok(run.front.length>3);
 for(const front of run.front){
   assert.equal(front.pass,true);
   assert.ok(run.results.every(other=>!M.dominates(other,front)));
 }
 const again=M.explore('ni',125,9172026),again2=M.explore('ni',125,9172026);
 assert.deepEqual(again.results.map(x=>x.asd1),again2.results.map(x=>x.asd1));
});
test('front-end script has parseable JS, references matching model script and one-axis mode',()=>{
 const html=fs.readFileSync(path.resolve(__dirname,'../fluxgate.html'),'utf8');
 assert.match(html,/src="fluxgate-model\.js"/);
 assert.match(html,/одна измерительная ось/);
 const scripts=[...html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)];
 assert.equal(scripts.length,1);
 for(const block of scripts)new vm.Script(block[1]);
 assert.ok(!html.includes('new THREE'));
 assert.match(html,/5000/);
});
