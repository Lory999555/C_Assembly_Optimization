%include "sseutils64.nasm"

section .data

subb equ 16

dim equ 4

section .bss

section .text

global rowDistance64AdcA

global rowDistance64AdcU

rowDistance64AdcA:

    push		rbp				; salva il Base Pointer
    mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
    pushaq						; salva i registri generali

    mov r11,[rbp+subb]
    ;vprintreg r11

    imul r8,dim                        ;prendo j
    ;vprintreg r8
    imul r8,r11                        ;j*sub
    ;vprintreg r8
    imul r8,r9                        ;j*k*subb
    ;vprintreg r8

    imul rcx,dim                        ;prendo i
    ;vprintreg rcx
    imul rcx,r11                        ;i*subb
    add rcx,r8                          ;j*k*subb+i*subb
    add rcx,rdi                         ;c[j*k*subb+i*subb]
    ;vprintreg rdx
    ;vprintreg rcx
    ;vxorps xmm5,xmm5                     ;distance=0 

    vmovaps ymm1,[rsi]                   ;uj_x[0]
    ;printregyps ymm1
    vmovaps ymm3,[rcx]
    ;printregyps ymm3
    vsubps ymm1,ymm3                   ;uj_x[0]-c[0*k*subb+i*subb]
    ;printregyps ymm1
    vmulps ymm1,ymm1                     ;(..)^2
    ;printregyps ymm1
    ;mov r11,[rbp+subb]                     ;subb
    sub r11,16                          ;subb-16

    mov r12,8                            ;z=8
    vxorps ymm2,ymm2
ciclo:
    cmp r12,r11
    jg resto

    vmovaps ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    vmovaps ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    jmp ciclo
resto:
    mov r11,[rbp+subb]                 ;subb
    sub r11,8                          ;subb-4
    ;vprintreg r11
ciclo2:
    cmp r12, r11
    jg fine
    vmovaps ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregyps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregyps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregyps ymm0
    vaddps ymm2,ymm0                     ;distance+=(..)^2
    ;printregyps ymm2
    add r12,8                           ;avanzo di indice
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

rowDistance64AdcU:

    push		rbp				; salva il Base Pointer
    mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
    pushaq						; salva i registri generali

    mov r11,[rbp+subb]
    ;vprintreg r11

    imul r8,dim                        ;prendo j
    ;vprintreg r8
    imul r8,r11                        ;j*sub
    ;vprintreg r8
    imul r8,r9                        ;j*k*subb
    ;vprintreg r8

    imul rcx,dim                        ;prendo i
    ;vprintreg rcx
    imul rcx,r11                        ;i*subb
    add rcx,r8                          ;j*k*subb+i*subb
    add rcx,rdi                         ;c[j*k*subb+i*subb]
    ;vprintreg rdx
    ;vprintreg rcx
    ;vxorps xmm5,xmm5                     ;distance=0 

    vmovups ymm1,[rsi]                   ;uj_x[0]
    ;printregyps ymm1
    vmovups ymm3,[rcx]
    ;printregyps ymm3
    vsubps ymm1,ymm3                   ;uj_x[0]-c[0*k*subb+i*subb]
    ;printregyps ymm1
    vmulps ymm1,ymm1                     ;(..)^2
    ;printregyps ymm1
    ;mov r11,[rbp+subb]                     ;subb
    sub r11,16                          ;subb-16

    mov r12,8                            ;z=8
    vxorps ymm2,ymm2
cicloU:
    cmp r12,r11
    jge restoU

    vmovups ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    vmovaps ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,8                           ;avanzo di indice

    jmp cicloU
restoU:
    mov r11,[rbp+subb]                 ;subb
    sub r11,8                          ;subb-4
    ;vprintreg r11
ciclo2U:
    cmp r12, r11
    jg resto2U
    vmovups ymm0,[rsi+4*r12]             ;uj_x[z]
    ;printregyps ymm0
    vsubps ymm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregyps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregyps ymm0
    vaddps ymm2,ymm0                     ;distance+=(..)^2
    ;printregyps ymm2
    add r12,8                           ;avanzo di indice
    jmp ciclo2U                       
resto2U:
    mov r11, [rbp+subb]
ciclo3U:
    cmp r12, r11
    je  fineU
    vmovss xmm0,[rsi+4*r12]             ;uj_x[z]
    ;printregps ymm0
    vsubss xmm0,[rcx+4*r12]              ;uj_x[z]- c[j*k*subb+i*subb+z]
    ;printregps ymm0
    vmulss xmm0,xmm0                     ;(..)^2
    ;printregps ymm0
    vaddss xmm2,xmm0                     ;distance+=(..)^2
    ;printregps ymm1
    add r12,1                           ;avanzo di indice
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