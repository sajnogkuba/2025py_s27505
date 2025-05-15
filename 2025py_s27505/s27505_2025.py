# Cel programu:
# Program generuje losową sekwencję DNA o długości podanej przez użytkownika, wstawia imię użytkownika w losowym miejscu,
# zapisuje wynik do pliku w formacie FASTA, a następnie wyświetla statystyki zawartości nukleotydów.

# Kontekst zastosowania:
# Program może być wykorzystywany w biologii molekularnej i bioinformatyce do testowania narzędzi analizujących dane DNA,
# generowania danych testowych, demonstracji formatu FASTA itp.

import random


# Funkcja do generowania losowej sekwencji DNA
def generate_dna_sequence(length):
    return ''.join(random.choices('ACGT', k=length))


# Funkcja do wstawienia imienia w losowe miejsce sekwencji
def insert_name(sequence, name):
    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name + sequence[pos:]


# Funkcja do liczenia statystyk nukleotydów
def calculate_stats(sequence):
    # Pomijamy litery niebędące A, C, G, T (czyli np. imię użytkownika)
    filtered = [nt for nt in sequence if nt in 'ACGT']
    total = len(filtered)
    stats = {nt: (filtered.count(nt) / total * 100) for nt in 'ACGT'}
    cg = stats['C'] + stats['G']
    at = stats['A'] + stats['T']
    cg_at_ratio = cg / at * 100 if at != 0 else 0
    return stats, cg_at_ratio


# Główna część programu
def main():
    # --- Pobieranie danych od użytkownika z walidacją ---
    # ORIGINAL:
    # length = int(input("Podaj długość sekwencji: "))
    # MODIFIED (dodano walidację wejścia – uodpornienie na błędne dane):
    while True:
        try:
            length = int(input("Podaj długość sekwencji: "))
            if length <= 0:
                raise ValueError
            break
        except ValueError:
            print("Wprowadź dodatnią liczbę całkowitą.")

    seq_id = input("Podaj ID sekwencji: ")
    description = input("Podaj opis sekwencji: ")
    name = input("Podaj imię: ")

    # --- Generowanie sekwencji i wstawianie imienia ---
    dna_seq = generate_dna_sequence(length)
    dna_with_name = insert_name(dna_seq, name)

    # --- Zapis do pliku FASTA ---
    filename = f"{seq_id}.fasta"
    with open(filename, 'w') as f:
        f.write(f">{seq_id} {description}\n")
        # ORIGINAL:
        # f.write(dna_with_name + '\n')
        # MODIFIED (podział sekwencji co 60 znaków zgodnie z konwencją FASTA):
        for i in range(0, len(dna_with_name), 60):
            f.write(dna_with_name[i:i + 60] + '\n')

    print(f"Sekwencja została zapisana do pliku {filename}")

    # --- Obliczanie i wyświetlanie statystyk ---
    stats, ratio = calculate_stats(dna_with_name)

    # ORIGINAL:
    # print("Statystyki sekwencji:")
    # for nt in 'ACGT':
    #     print(f"{nt}: {stats[nt]:.1f}%")
    # print(f"%CG: {ratio:.1f}")
    # MODIFIED (dodano nagłówek i nazwę stosunku):
    print("\n--- Statystyki sekwencji ---")
    for nt in 'ACGT':
        print(f"{nt}: {stats[nt]:.1f}%")
    print(f"Stosunek CG/AT: {ratio:.1f}%")


# Uruchomienie programu
if __name__ == "__main__":
    main()
