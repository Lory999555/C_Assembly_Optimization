%include "sseutils64.nasm"

section .data
i equ 8
j equ 12
n equ 16


indice equ 20


index equ 24

section .bss

section .text

global mapping32

mapping32:
    push	ebp							; salva il Base Pointer
    mov		ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
    push	ebx							; salva i registri da preservare
    push	esi
    push	edi

    mov ebx,[ebp+i]
    mov ecx,[ebp+j]

    cmp ebx,ecx
    je zero 

    mov eax,[ebp+n]
    mov esi,[ebp+index]    ;n*(n-1)/2 in esi

    cmp ebx,ecx
    jg iMag

    sub eax,ebx ;n-i
    mov edi,eax ;n-i
    sub edi,1   ;n-i-1
    imul eax,edi ;(n-i-1)*(n-i)
    mov edi,2
    xor edx,edx
    div edi      ;(n-i-1)*(n-j)/2

    sub esi,eax   ;n*(n-1)/2-(n-i-1)*(n-j)/2

    add esi,ecx ;n*(n-1)/2-(n-i-1)*(n-j)/2+j
    sub esi,ebx ;n*(n-1)/2-(n-i-1)*(n-j)/2+j-i
    sub esi,1   ;n*(n-1)/2-(n-i-1)*(n-j)/2+j-i-1
    jmp fine
iMag:
    sub eax,ecx ;n-j
    mov edi,eax ;n-j
    sub edi,1   ;n-j-1
    imul eax,edi ;(n-j-1)*(n-j)
    mov edi,2
    xor edx,edx
    div edi      ;(n-j-1)*(n-j)/2
    
    sub esi,eax   ;n*(n-1)/2-(n-j-1)*(n-j)/2

    add esi,ebx ;n*(n-1)/2-(n-i-1)*(n-j)/2+i
    sub esi,ecx ;n*(n-1)/2-(n-i-1)*(n-j)/2+i-j
    sub esi,1   ;n*(n-1)/2-(n-i-1)*(n-j)/2+i-j-1
    jmp fine
zero:
    mov esi,-1
fine:
    mov ebx,[ebp+indice]
    mov [ebx],esi
    pop	edi									; ripristina i registri da preservare
    pop	esi
    pop	ebx
    mov	esp, ebp							; ripristina lo Stack Pointer
    pop	ebp									; ripristina il Base Pointer
    ret	