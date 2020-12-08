{
  let col = ref 1

  let program = ref ""

  let string_buffer = Buffer.create 1024

  let process_label = fun s ->
    let size = String.length s in
    col := !col + size;
    program := !program ^ s

  let process_instruction = fun s ->
    let size = String.length s in
    for i = !col to 16 do
      incr col;
      program := !program ^ " "
    done;
    col := !col + size;
    program := !program ^ s

  let process_operand = fun s ->
    for i = !col to 26 do
      incr col;
      program := !program ^ " "
    done;
    let size = String.length s in
    col := !col + size;
    program := !program ^ s

  let process_comment =
    fun s ->
      let size = String.length s in
      if !col = 1 then
        begin
          col := !col + size;
          program := !program ^ s
        end
      else
        begin
          for i = !col to 80 do
            incr col;
            program := !program ^ " "
          done;
          col := !col + size;
          program := !program ^ s
        end

  let process_identifier =
    let mnemonic = Hashtbl.create 67 in
    (* NASM *)
    List.iter (fun s -> Hashtbl.add mnemonic s ()) [
      "global"; "extern"; "section"; "segment"; "org"; "db"; "dw"; "dd"; "dq"; "equ"; "resb"; "resw"; "resd"; "resq"; "je"; "jz"; "jne"; "jnz"; "jg"; "jnle"; "jge"; "jnl"; "jl"; "jnge"; "jle"; "jng"; "je"; "jz"; "jne"; "jnz"; "ja"; "jnbe"; "jae"; "jnb"; "jb"; "jnae"; "jbe"; "jna"; "jxcz"; "jc"; "jnc"; "jo"; "jno"; "jp"; "jpe"; "jnp"; "jpo"; "js"; "jns";
    ];
    (* x86_64 instructions *)
    List.iter (fun s -> Hashtbl.add mnemonic s ()) [
      "aaa"; "aad"; "aam"; "aas"; "adc"; "adcx"; "add"; "addpd"; "addps"; "addsd"; "addss"; "addsubpd"; "addsubps"; "adox"; "aesdec"; "aesdeclast"; "aesenc"; "aesenclast"; "aesimc"; "aeskeygenassist"; "and"; "andn"; "andnpd"; "andnps"; "andpd"; "andps"; "arpl"; "bextr"; "blendpd"; "blendps"; "blendvpd"; "blendvps"; "blsi"; "blsmsk"; "blsr"; "bndcl"; "bndcn"; "bndcu"; "bndldx"; "bndmk"; "bndmov"; "bndstx"; "bound"; "bsf"; "bsr"; "bswap"; "bt"; "btc"; "btr"; "bts"; "bzhi"; "call"; "cbw"; "cdq"; "cdqe"; "clac"; "clc"; "cld"; "cldemote"; "clflush"; "clflushopt"; "cli"; "clts"; "clwb"; "cmc"; "cmovcc"; "cmp"; "cmppd"; "cmpps"; "cmps"; "cmpsb"; "cmpsd"; "cmpsd"; "cmpsq"; "cmpss"; "cmpsw"; "cmpxchg"; "cmpxchg16b"; "cmpxchg8b"; "comisd"; "comiss"; "cpuid"; "cqo"; "crc32"; "cvtdq2pd"; "cvtdq2ps"; "cvtpd2dq"; "cvtpd2pi"; "cvtpd2ps"; "cvtpi2pd"; "cvtpi2ps"; "cvtps2dq"; "cvtps2pd"; "cvtps2pi"; "cvtsd2si"; "cvtsd2ss"; "cvtsi2sd"; "cvtsi2ss"; "cvtss2sd"; "cvtss2si"; "cvttpd2dq"; "cvttpd2pi"; "cvttps2dq"; "cvttps2pi"; "cvttsd2si"; "cvttss2si"; "cwd"; "cwde"; "daa"; "das"; "dec"; "div"; "divpd"; "divps"; "divsd"; "divss"; "dppd"; "dpps"; "emms"; "enter"; "extractps"; "f2xm1"; "fabs"; "fadd"; "faddp"; "fbld"; "fbstp"; "fchs"; "fclex"; "fcmovcc"; "fcom"; "fcomi"; "fcomip"; "fcomp"; "fcompp"; "fcos"; "fdecstp"; "fdiv"; "fdivp"; "fdivr"; "fdivrp"; "ffree"; "fiadd"; "ficom"; "ficomp"; "fidiv"; "fidivr"; "fild"; "fimul"; "fincstp"; "finit"; "fist"; "fistp"; "fisttp"; "fisub"; "fisubr"; "fld"; "fld1"; "fldcw"; "fldenv"; "fldl2e"; "fldl2t"; "fldlg2"; "fldln2"; "fldpi"; "fldz"; "fmul"; "fmulp"; "fnclex"; "fninit"; "fnop"; "fnsave"; "fnstcw"; "fnstenv"; "fnstsw"; "fpatan"; "fprem"; "fprem1"; "fptan"; "frndint"; "frstor"; "fsave"; "fscale"; "fsin"; "fsincos"; "fsqrt"; "fst"; "fstcw"; "fstenv"; "fstp"; "fstsw"; "fsub"; "fsubp"; "fsubr"; "fsubrp"; "ftst"; "fucom"; "fucomi"; "fucomip"; "fucomp"; "fucompp"; "fwait"; "fxam"; "fxch"; "fxrstor"; "fxsave"; "fxtract"; "fyl2x"; "fyl2xp1"; "gf2p8affineinvqb"; "gf2p8affineqb"; "gf2p8mulb"; "haddpd"; "haddps"; "hlt"; "hsubpd"; "hsubps"; "idiv"; "imul"; "in"; "inc"; "ins"; "insb"; "insd"; "insertps"; "insw"; "int"; "int1"; "int3"; "into"; "invd"; "invlpg"; "invpcid"; "iret"; "iretd"; "jmp"; "jcc"; "kaddb"; "kaddd"; "kaddq"; "kaddw"; "kandb"; "kandd"; "kandnb"; "kandnd"; "kandnq"; "kandnw"; "kandq"; "kandw"; "kmovb"; "kmovd"; "kmovq"; "kmovw"; "knotb"; "knotd"; "knotq"; "knotw"; "korb"; "kord"; "korq"; "kortestb"; "kortestd"; "kortestq"; "kortestw"; "korw"; "kshiftlb"; "kshiftld"; "kshiftlq"; "kshiftlw"; "kshiftrb"; "kshiftrd"; "kshiftrq"; "kshiftrw"; "ktestb"; "ktestd"; "ktestq"; "ktestw"; "kunpckbw"; "kunpckdq"; "kunpckwd"; "kxnorb"; "kxnord"; "kxnorq"; "kxnorw"; "kxorb"; "kxord"; "kxorq"; "kxorw"; "lahf"; "lar"; "lddqu"; "ldmxcsr"; "lds"; "lea"; "leave"; "les"; "lfence"; "lfs"; "lgdt"; "lgs"; "lidt"; "lldt"; "lmsw"; "lock"; "lods"; "lodsb"; "lodsd"; "lodsq"; "lodsw"; "loop"; "loopcc"; "lsl"; "lss"; "ltr"; "lzcnt"; "maskmovdqu"; "maskmovq"; "maxpd"; "maxps"; "maxsd"; "maxss"; "mfence"; "minpd"; "minps"; "minsd"; "minss"; "monitor"; "mov"; "mov"; "mov"; "movapd"; "movaps"; "movbe"; "movd"; "movddup"; "movdir64b"; "movdiri"; "movdq2q"; "movdqa"; "movdqu"; "movhlps"; "movhpd"; "movhps"; "movlhps"; "movlpd"; "movlps"; "movmskpd"; "movmskps"; "movntdq"; "movntdqa"; "movnti"; "movntpd"; "movntps"; "movntq"; "movq"; "movq"; "movq2dq"; "movs"; "movsb"; "movsd"; "movsd"; "movshdup"; "movsldup"; "movsq"; "movss"; "movsw"; "movsx"; "movsxd"; "movupd"; "movups"; "movzx"; "mpsadbw"; "mul"; "mulpd"; "mulps"; "mulsd"; "mulss"; "mulx"; "mwait"; "neg"; "nop"; "not"; "or"; "orpd"; "orps"; "out"; "outs"; "outsb"; "outsd"; "outsw"; "pabsb"; "pabsd"; "pabsq"; "pabsw"; "packssdw"; "packsswb"; "packusdw"; "packuswb"; "paddb"; "paddd"; "paddq"; "paddsb"; "paddsw"; "paddusb"; "paddusw"; "paddw"; "palignr"; "pand"; "pandn"; "pause"; "pavgb"; "pavgw"; "pblendvb"; "pblendw"; "pclmulqdq"; "pcmpeqb"; "pcmpeqd"; "pcmpeqq"; "pcmpeqw"; "pcmpestri"; "pcmpestrm"; "pcmpgtb"; "pcmpgtd"; "pcmpgtq"; "pcmpgtw"; "pcmpistri"; "pcmpistrm"; "pdep"; "pext"; "pextrb"; "pextrd"; "pextrq"; "pextrw"; "phaddd"; "phaddsw"; "phaddw"; "phminposuw"; "phsubd"; "phsubsw"; "phsubw"; "pinsrb"; "pinsrd"; "pinsrq"; "pinsrw"; "pmaddubsw"; "pmaddwd"; "pmaxsb"; "pmaxsd"; "pmaxsq"; "pmaxsw"; "pmaxub"; "pmaxud"; "pmaxuq"; "pmaxuw"; "pminsb"; "pminsd"; "pminsq"; "pminsw"; "pminub"; "pminud"; "pminuq"; "pminuw"; "pmovmskb"; "pmovsx"; "pmovzx"; "pmuldq"; "pmulhrsw"; "pmulhuw"; "pmulhw"; "pmulld"; "pmullq"; "pmullw"; "pmuludq"; "pop"; "popa"; "popad"; "popcnt"; "popf"; "popfd"; "popfq"; "por"; "prefetchw"; "prefetchh"; "psadbw"; "pshufb"; "pshufd"; "pshufhw"; "pshuflw"; "pshufw"; "psignb"; "psignd"; "psignw"; "pslld"; "pslldq"; "psllq"; "psllw"; "psrad"; "psraq"; "psraw"; "psrld"; "psrldq"; "psrlq"; "psrlw"; "psubb"; "psubd"; "psubq"; "psubsb"; "psubsw"; "psubusb"; "psubusw"; "psubw"; "ptest"; "ptwrite"; "punpckhbw"; "punpckhdq"; "punpckhqdq"; "punpckhwd"; "punpcklbw"; "punpckldq"; "punpcklqdq"; "punpcklwd"; "push"; "pusha"; "pushad"; "pushf"; "pushfd"; "pushfq"; "pxor"; "rcl"; "rcpps"; "rcpss"; "rcr"; "rdfsbase"; "rdgsbase"; "rdmsr"; "rdpid"; "rdpkru"; "rdpmc"; "rdrand"; "rdseed"; "rdtsc"; "rdtscp"; "rep"; "repe"; "repne"; "repnz"; "repz"; "ret"; "rol"; "ror"; "rorx"; "roundpd"; "roundps"; "roundsd"; "roundss"; "rsm"; "rsqrtps"; "rsqrtss"; "sahf"; "sal"; "sar"; "sarx"; "sbb"; "scas"; "scasb"; "scasd"; "scasw"; "setcc"; "sfence"; "sgdt"; "sha1msg1"; "sha1msg2"; "sha1nexte"; "sha1rnds4"; "sha256msg1"; "sha256msg2"; "sha256rnds2"; "shl"; "shld"; "shlx"; "shr"; "shrd"; "shrx"; "shufpd"; "shufps"; "sidt"; "sldt"; "smsw"; "sqrtpd"; "sqrtps"; "sqrtsd"; "sqrtss"; "stac"; "stc"; "std"; "sti"; "stmxcsr"; "stos"; "stosb"; "stosd"; "stosq"; "stosw"; "str"; "sub"; "subpd"; "subps"; "subsd"; "subss"; "swapgs"; "syscall"; "sysenter"; "sysexit"; "sysret"; "test"; "tpause"; "tzcnt"; "ucomisd"; "ucomiss"; "ud"; "umonitor"; "umwait"; "unpckhpd"; "unpckhps"; "unpcklpd"; "unpcklps"; "valignd"; "valignq"; "vblendmpd"; "vblendmps"; "vbroadcast"; "vcompresspd"; "vcompressps"; "vcvtpd2qq"; "vcvtpd2udq"; "vcvtpd2uqq"; "vcvtph2ps"; "vcvtps2ph"; "vcvtps2qq"; "vcvtps2udq"; "vcvtps2uqq"; "vcvtqq2pd"; "vcvtqq2ps"; "vcvtsd2usi"; "vcvtss2usi"; "vcvttpd2qq"; "vcvttpd2udq"; "vcvttpd2uqq"; "vcvttps2qq"; "vcvttps2udq"; "vcvttps2uqq"; "vcvttsd2usi"; "vcvttss2usi"; "vcvtudq2pd"; "vcvtudq2ps"; "vcvtuqq2pd"; "vcvtuqq2ps"; "vcvtusi2sd"; "vcvtusi2ss"; "vdbpsadbw"; "verr"; "verw"; "vexpandpd"; "vexpandps"; "vextractf128"; "vextractf32x4"; "vextractf32x8"; "vextractf64x2"; "vextractf64x4"; "vextracti128"; "vextracti32x4"; "vextracti32x8"; "vextracti64x2"; "vextracti64x4"; "vfixupimmpd"; "vfixupimmps"; "vfixupimmsd"; "vfixupimmss"; "vfmadd132pd"; "vfmadd132ps"; "vfmadd132sd"; "vfmadd132ss"; "vfmadd213pd"; "vfmadd213ps"; "vfmadd213sd"; "vfmadd213ss"; "vfmadd231pd"; "vfmadd231ps"; "vfmadd231sd"; "vfmadd231ss"; "vfmaddsub132pd"; "vfmaddsub132ps"; "vfmaddsub213pd"; "vfmaddsub213ps"; "vfmaddsub231pd"; "vfmaddsub231ps"; "vfmsub132pd"; "vfmsub132ps"; "vfmsub132sd"; "vfmsub132ss"; "vfmsub213pd"; "vfmsub213ps"; "vfmsub213sd"; "vfmsub213ss"; "vfmsub231pd"; "vfmsub231ps"; "vfmsub231sd"; "vfmsub231ss"; "vfmsubadd132pd"; "vfmsubadd132ps"; "vfmsubadd213pd"; "vfmsubadd213ps"; "vfmsubadd231pd"; "vfmsubadd231ps"; "vfnmadd132pd"; "vfnmadd132ps"; "vfnmadd132sd"; "vfnmadd132ss"; "vfnmadd213pd"; "vfnmadd213ps"; "vfnmadd213sd"; "vfnmadd213ss"; "vfnmadd231pd"; "vfnmadd231ps"; "vfnmadd231sd"; "vfnmadd231ss"; "vfnmsub132pd"; "vfnmsub132ps"; "vfnmsub132sd"; "vfnmsub132ss"; "vfnmsub213pd"; "vfnmsub213ps"; "vfnmsub213sd"; "vfnmsub213ss"; "vfnmsub231pd"; "vfnmsub231ps"; "vfnmsub231sd"; "vfnmsub231ss"; "vfpclasspd"; "vfpclassps"; "vfpclasssd"; "vfpclassss"; "vgatherdpd"; "vgatherdpd"; "vgatherdps"; "vgatherdps"; "vgatherqpd"; "vgatherqpd"; "vgatherqps"; "vgatherqps"; "vgetexppd"; "vgetexpps"; "vgetexpsd"; "vgetexpss"; "vgetmantpd"; "vgetmantps"; "vgetmantsd"; "vgetmantss"; "vinsertf128"; "vinsertf32x4"; "vinsertf32x8"; "vinsertf64x2"; "vinsertf64x4"; "vinserti128"; "vinserti32x4"; "vinserti32x8"; "vinserti64x2"; "vinserti64x4"; "vmaskmov"; "vmovdqa32"; "vmovdqa64"; "vmovdqu16"; "vmovdqu32"; "vmovdqu64"; "vmovdqu8"; "vpblendd"; "vpblendmb"; "vpblendmd"; "vpblendmq"; "vpblendmw"; "vpbroadcast"; "vpbroadcastb"; "vpbroadcastd"; "vpbroadcastm"; "vpbroadcastq"; "vpbroadcastw"; "vpcmpb"; "vpcmpd"; "vpcmpq"; "vpcmpub"; "vpcmpud"; "vpcmpuq"; "vpcmpuw"; "vpcmpw"; "vpcompressd"; "vpcompressq"; "vpconflictd"; "vpconflictq"; "vperm2f128"; "vperm2i128"; "vpermb"; "vpermd"; "vpermi2b"; "vpermi2d"; "vpermi2pd"; "vpermi2ps"; "vpermi2q"; "vpermi2w"; "vpermilpd"; "vpermilps"; "vpermpd"; "vpermps"; "vpermq"; "vpermt2b"; "vpermt2d"; "vpermt2pd"; "vpermt2ps"; "vpermt2q"; "vpermt2w"; "vpermw"; "vpexpandd"; "vpexpandq"; "vpgatherdd"; "vpgatherdd"; "vpgatherdq"; "vpgatherdq"; "vpgatherqd"; "vpgatherqd"; "vpgatherqq"; "vpgatherqq"; "vplzcntd"; "vplzcntq"; "vpmadd52huq"; "vpmadd52luq"; "vpmaskmov"; "vpmovb2m"; "vpmovd2m"; "vpmovdb"; "vpmovdw"; "vpmovm2b"; "vpmovm2d"; "vpmovm2q"; "vpmovm2w"; "vpmovq2m"; "vpmovqb"; "vpmovqd"; "vpmovqw"; "vpmovsdb"; "vpmovsdw"; "vpmovsqb"; "vpmovsqd"; "vpmovsqw"; "vpmovswb"; "vpmovusdb"; "vpmovusdw"; "vpmovusqb"; "vpmovusqd"; "vpmovusqw"; "vpmovuswb"; "vpmovw2m"; "vpmovwb"; "vpmultishiftqb"; "vprold"; "vprolq"; "vprolvd"; "vprolvq"; "vprord"; "vprorq"; "vprorvd"; "vprorvq"; "vpscatterdd"; "vpscatterdq"; "vpscatterqd"; "vpscatterqq"; "vpsllvd"; "vpsllvq"; "vpsllvw"; "vpsravd"; "vpsravq"; "vpsravw"; "vpsrlvd"; "vpsrlvq"; "vpsrlvw"; "vpternlogd"; "vpternlogq"; "vptestmb"; "vptestmd"; "vptestmq"; "vptestmw"; "vptestnmb"; "vptestnmd"; "vptestnmq"; "vptestnmw"; "vrangepd"; "vrangeps"; "vrangesd"; "vrangess"; "vrcp14pd"; "vrcp14ps"; "vrcp14sd"; "vrcp14ss"; "vreducepd"; "vreduceps"; "vreducesd"; "vreducess"; "vrndscalepd"; "vrndscaleps"; "vrndscalesd"; "vrndscaless"; "vrsqrt14pd"; "vrsqrt14ps"; "vrsqrt14sd"; "vrsqrt14ss"; "vscalefpd"; "vscalefps"; "vscalefsd"; "vscalefss"; "vscatterdpd"; "vscatterdps"; "vscatterqpd"; "vscatterqps"; "vshuff32x4"; "vshuff64x2"; "vshufi32x4"; "vshufi64x2"; "vtestpd"; "vtestps"; "vzeroall"; "vzeroupper"; "wait"; "wbinvd"; "wrfsbase"; "wrgsbase"; "wrmsr"; "wrpkru"; "xabort"; "xacquire"; "xadd"; "xbegin"; "xchg"; "xend"; "xgetbv"; "xlat"; "xlatb"; "xor"; "xorpd"; "xorps"; "xrelease"; "xrstor"; "xrstors"; "xsave"; "xsavec"; "xsaveopt"; "xsaves"; "xsetbv"; "xtest";
    ];
    fun s ->
      try
        Hashtbl.find mnemonic s;
        (* instruction *)
        process_instruction s
      with Not_found ->
        (* operand *)
        process_operand s
}

let space = [' ' '\t']+

let letter = ['a'-'z' 'A'-'Z']

let digit = ['0'-'9']

let number = "0x"? digit+

let name =  '_'* letter+ '_'* | '_' digit+ '_'*

let ident = '.'? name+ | (name '_'+ name)+

let label = ident ':'

let comment = ";" [^ '\n']*

rule format = parse
| space             { format lexbuf }
| '\n'              { program := !program ^ "\n"; col := 1; format lexbuf }
| '"'               { process_operand (quotes lexbuf); format lexbuf }
| '['               { process_operand (brackets lexbuf); format lexbuf }
| ','               { program := !program ^ ", "; col := !col + 2; format lexbuf }
| label as lbl      { process_label lbl; format lexbuf }
| ident as id       { process_identifier id; format lexbuf }
| number as n       { process_operand n; format lexbuf }
| comment as cmt    { process_comment cmt; format lexbuf }
| _ as c            { program := !program ^ (String.make 1 c); incr col; format lexbuf }
| eof               { !program }

and format_no_comments = parse
| space             { format_no_comments lexbuf }
| '\n'              { program := !program ^ "\n"; col := 1; format_no_comments lexbuf }
| '"'               { process_operand (quotes lexbuf); format_no_comments lexbuf }
| '['               { process_operand (brackets lexbuf); format_no_comments lexbuf }
| ','               { program := !program ^ ", "; col := !col + 2; format_no_comments lexbuf }
| label as lbl      { process_label lbl; format_no_comments lexbuf }
| ident as id       { process_identifier id; format_no_comments lexbuf }
| number as n       { process_operand n; format_no_comments lexbuf }
| comment           { format_no_comments lexbuf }
| _ as c            { program := !program ^ (String.make 1 c); incr col; format_no_comments lexbuf }
| eof               { !program }

and quotes = parse
  | '"'
    {
      let s = "\"" ^ Buffer.contents string_buffer ^ "\"" in
      Buffer.reset string_buffer;
      s
    }
  | _ as c
    {
      Buffer.add_char string_buffer c;
      quotes lexbuf
    }

and brackets = parse
  | ']'
    {
      let s = "[" ^ Buffer.contents string_buffer ^ "]" in
      Buffer.reset string_buffer;
      s
    }
  | ident as id
    {
      String.iter (fun c -> Buffer.add_char string_buffer c) id;
      brackets lexbuf
    }
  | _ as c
    {
      Buffer.add_char string_buffer c;
      brackets lexbuf
    }

{
  let main () =
    let ifile = ref "" in
    let ofile = ref "" in
    let inplace = ref false in
    let comments = ref false in
    let help = ref false in

    let option = [
      "-i", Arg.Set inplace, "";
      "-o", Arg.String (fun s -> ofile := s), "";
      "--no-comments", Arg.Set comments, "";
      "-h", Arg.Set help, "";
      "--help", Arg.Set help, "";
    ]
    in

    let usage = "usage: fmtnasm [options] [filename]" in
    let info =
      usage ^ "\n\n" ^
      "-i\t\t\tformat in place\n" ^
      "-o <file>\t\tselect output file\n\n" ^
      "--no-comments\t\tremove comments flag\n\n" ^
      "-h\t\t\tprint help\n" ^
      "--help\t\t\tprint help\n"
    in

    let set_input s = ifile := s in
    Arg.parse option set_input usage;

    if !help then
      begin
        Printf.printf "%s" info;
        exit 0
      end;

    let lexbuf =
      if !ifile <> "" then
        let c = open_in !ifile in
        Lexing.from_channel c
      else
        Lexing.from_channel stdin
    in

    if !inplace && !ifile <> "" then
      ofile := !ifile;

    let program =
      if !comments then
        format_no_comments lexbuf
      else
        format lexbuf
    in

    let print =
      if !ofile <> "" then
        let c = open_out !ofile in
        Printf.fprintf c "%s"
      else
        Printf.printf "%s"
    in

    print program

  let _ = Printexc.print main ()
}
