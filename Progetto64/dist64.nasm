%include "sseutils64.nasm"

section .data

dim equ 4

section .bss

section .text

global dist64A

dist64A:
    start

    vmovaps ymm1,[rdi] ;x[0]
    vsubps ymm1,[rsi]  ;x[0]-y[0]
    vmulps ymm1,ymm1     ;(..)^2
    ;printregyps ymm1
    mov r12,rcx         ;d
    ;mov ecx,dim     ;4
    sub r12,16     ;d-16

    vxorps ymm2, ymm2
    mov r10,8       ;i=4
ciclo:
    cmp r10,r12     ;(j>=d-16)?
    jg resto
    vmovaps ymm0,[rdi+4*r10]
    vmovaps ymm7,[rdi+4*r10+32] ;x[i]
    vsubps ymm0,[rsi+4*r10]  ;x[i]-y[i]
    vsubps ymm7,[rsi+4*r10+32]
    vmulps ymm0,ymm0         ;(..)^2
    vmulps ymm7,ymm7
    ;printregyps ymm0
    vaddps ymm1,ymm0         ;distance+=(..)^2
    vaddps ymm1,ymm7
    ;printregyps ymm1
    add r10,8               ;avanzo di indice

    ;printregyps ymm0
    ;printregyps ymm1
    add r10,8

    jmp ciclo
resto:
    mov r12, rcx
    sub r12, 8
ciclo2:
    cmp r10, r12
    jg  fine
    vmovaps ymm0,[rdi+4*r10] ;sommo gli ultimi elementi rimanenti
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
    vmovss  [rdx],xmm1
    stop

global dist64U

dist64U:
    start

    vmovups ymm1,[rdi] ;x[0]
    vsubps ymm1,[rsi]  ;x[0]-y[0]
    vmulps ymm1,ymm1     ;(..)^2
    ;printregyps ymm1
    mov r12,rcx         ;d
    ;mov ecx,dim     ;4
    sub r12,16     ;d-16

    vxorps ymm2, ymm2
    mov r10,8       ;i=4
cicloU:
    cmp r10,r12     ;(j>=d-16)?
    jg restoU
    vmovups ymm0,[rdi+4*r10] ;x[i]
    vmovups ymm7,[rdi+4*r10+32]
    vsubps ymm0,[rsi+4*r10]  ;x[i]-y[i]
    vsubps ymm7,[rsi+4*r10+32]
    vmulps ymm0,ymm0         ;(..)^2
    vmulps ymm7,ymm7
    ;printregyps ymm0
    vaddps ymm1,ymm0         ;distance+=(..)^2
    vaddps ymm1,ymm7
    ;printregyps ymm1
    add r10,8               ;avanzo di indice

    ;printregyps ymm0
    ;printregyps ymm1
    add r10,8

    jmp cicloU
restoU:
    mov r12, rcx
    sub r12, 8
ciclo2U:
    cmp r10, r12
    jg  resto2U
    vmovups ymm0,[rdi+4*r10] ;sommo gli ultimi elementi rimanenti
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
    vmovss xmm7,[rdi+4*r10] ;sommo gli ultimi elementi rimanenti
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
    vmovss  [rdx],xmm1
    stop

