%include "sseutils64.nasm"

section .data

dis     equ     16  ;dis temporanea
dim1    equ     4

section .bss
section .text

global cent_XA

cent_XA:
    start
    mov     r14, [rbp+dis]  ;dis

    mov     r11, 1
    ;sub     rdx, 4      ;k-4

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    
    vxorps   xmm3, xmm3
    vmovss   xmm3, [r14]         ;dis
foriA: 
    cmp     r11, rdx            ; i> k-4
    jg      fine_A

    mov     r13, rcx            ;d
    imul    r13,r11, dim1      ;4*i*d
    vprintreg r13
    vprintreg rdi
    add     r13,rdi            ;cent[4*i*d]
    vprintreg r13

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti1A
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i

avanti1A:
    add     r11, 1

    vprintreg r13
    vprintreg rdi
    add     r13,rdi            ;cent[4*i*d]
    vprintreg rdi

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti2A
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i
avanti2A:
    add     r11, 1

    add     rdi,r13             ;cent[4*i*d]
    vprintreg rdi

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti3A
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i

avanti3A:
    add     r11, 1

    add     rdi,r13             ;cent[4*i*d]
    vprintreg rdi

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64A
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jl      avanti4A
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i
avanti4A:
    add     r11, 1         
    jmp     foriA
fine_A:
    stop


global cent_XU

cent_XU:
    start
    mov     r14, [rbp+dis]  ;dis

    mov     r11, 0
    sub     rdx, 4      ;k-4

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    
    vxorps   xmm3, xmm3
    vmovss   xmm3, [r14]         ;dis
foriU: 
    cmp     r11, rdx            ; i> k-4
    jg      fine_U
    vprintreg r11

    mov     r13, rcx            ;d
    imul    r13,r11, dim1      ;4*i*d
    add     rdi,r13             ;cent[4*i*d]

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti1U
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i

avanti1U:
    add     r11, 1

    mov     r13, rcx            ;d
    imul    r13, r11            ;i*d
    imul    r13, dim1           ;4*i*d
    add     rdi,r13             ;cent[4*i*d]

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti2U
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i
avanti2U:
    add     r11, 1

    mov     r13, rcx            ;d
    imul    r13, r11            ;i*d
    imul    r13, dim1           ;4*i*d
    add     rdi,r13             ;cent[4*i*d]

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jb      avanti3U
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i

avanti3U:
    add     r11, 1

    mov     r13, rcx            ;d
    imul    r13, r11            ;i*d
    imul    r13, dim1           ;4*i*d
    add     rdi,r13             ;cent[4*i*d]

    push    rdi
    push    r13
    push    r12
    push    r10
    call    dist64U
    pop     r10
    pop     r12
    pop     r13
    pop     rdi

    ;vxorps   xmm4, xmm4      
    vmovss   xmm4, [r14]     ;dis
    vcomiss  xmm3, xmm4      ;tmp<dis?
    jl      avanti4U
    vmovss   xmm3, xmm4      ;dis=tmp
    mov     [r9], r11       ;park=i
avanti4U:
    add     r11, 1         
    jmp     foriU
fine_U:
    stop

section .data

dim equ 4

section .bss

section .text

dist64A:
    start
    
    vmovaps ymm1,[rdi] ;x[0]
    vsubps ymm1,[rsi]  ;x[0]-y[0]
    vmulps ymm1,ymm1     ;(..)^2
    ;printregyps ymm1
    mov r12,rcx         ;d
    ;mov rcx,dim     ;4
    sub r12,16     ;d-16

    vxorps ymm2, ymm2
    mov r10,8       ;i=4
ciclo:
    cmp r10,r12     ;(j>=d-16)?
    jg resto
    vmovaps ymm0, [rdi+4*r10] ;x[i]
    vsubps ymm0,[rsi+4*r10]  ;x[i]-y[i]
    vmulps ymm0,ymm0         ;(..)^2
    ;printregyps ymm0
    vaddps ymm1,ymm0         ;distance+=(..)^2
    ;printregyps ymm1
    add r10,8               ;avanzo di indice

    vmovaps ymm0, [rdi+4*r10]
    vsubps ymm0,[rsi+4*r10]
    vmulps ymm0,ymm0
    ;printregyps ymm0
    vaddps ymm1,ymm0
    ;printregyps ymm1
    add r10,8

    jmp ciclo
resto:
    mov r12, rcx
    sub r12, 8
ciclo2:
    cmp r10, r12
    jg  fine
    vmovaps ymm0, [rdi+4*r10] ;sommo gli ultimi elementi rimanenti
    vsubps ymm0,[rsi+4*r10]
    vmulps ymm0,ymm0
    vaddps ymm2, ymm0
    add r10, 8
    jmp ciclo2

fine:
    vaddps ymm1,ymm2        ;merge di tutte le somme
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vperm2f128 ymm5,ymm1,ymm1,1
    vaddss xmm1,xmm5
    ;printregyps ymm1
    vmovss  [r14],xmm1
    stop

dist64U:
    start

    vmovups ymm1,[rdi] ;x[0]
    vsubps ymm1,[rsi]  ;x[0]-y[0]
    vmulps ymm1,ymm1     ;(..)^2
    ;printregyps ymm1
    mov r12,rcx         ;d
    ;mov rcx,dim     ;4
    sub r12,16     ;d-16

    vxorps ymm2, ymm2
    mov r10,8       ;i=4
cicloU:
    cmp r10,r12     ;(j>=d-16)?
    jg restoU
    vmovups ymm0, [rdi+4*r10] ;x[i]
    vsubps ymm0,[rsi+4*r10]  ;x[i]-y[i]
    vmulps ymm0,ymm0         ;(..)^2
    ;printregyps ymm0
    vaddps ymm1,ymm0         ;distance+=(..)^2
    ;printregyps ymm1
    add r10,8               ;avanzo di indice

    vmovups ymm0, [rdi+4*r10]
    vsubps ymm0,[rsi+4*r10]
    vmulps ymm0,ymm0
    ;printregyps ymm0
    vaddps ymm1,ymm0
    ;printregyps ymm1
    add r10,8

    jmp cicloU
restoU:
    mov r12, rcx
    sub r12, 8
ciclo2U:
    cmp r10, r12
    jg  resto2U
    vmovups ymm0, [rdi+4*r10] ;sommo gli ultimi elementi rimanenti
    vsubps ymm0,[rsi+4*r10]
    vmulps ymm0,ymm0
    vaddps ymm2, ymm0
    ;printregyps ymm2
    add r10, 8
    jmp ciclo2U

resto2U:
    mov r12, rcx
ciclo3U:
    cmp r10, r12
    je  fineU
    vmovss xmm7, [rdi+4*r10] ;sommo gli ultimi elementi rimanenti
    vsubss xmm7,[rsi+4*r10]
    vmulss xmm7, xmm7
    ;printregyps ymm7
    vaddps ymm2, ymm7
    ;printregyps ymm2
    add r10, 1
    jmp ciclo3U
fineU:
    vaddps ymm1,ymm2        ;merge di tutte le somme
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vperm2f128 ymm5,ymm1,ymm1,1
    vaddss xmm1,xmm5
    ;printregyps ymm1
    vmovss  [r14],xmm1
    stop